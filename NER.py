import pandas as pd
import spacy
from spacy.pipeline import Sentencizer
from scispacy.abbreviation import AbbreviationDetector
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md

# Load the models
print("Loading models...")
nlp_bi = en_ner_bionlp13cg_md.load()
nlp_bc = en_ner_bc5cdr_md.load()
nlp_custom = spacy.load('./model-best')  # Load your custom model here
nlp_edc_custom = spacy.load('./EDC_target_model/output/EDC_model-best')
print("Models loaded successfully.")

# Add AbbreviationDetector to each model's pipeline
print("Adding AbbreviationDetector...")
for nlp in [nlp_bi, nlp_bc, nlp_custom]:
    nlp.add_pipe("abbreviation_detector")
print("AbbreviationDetector added successfully.")

# Add a sentencizer to the custom model if it does not have sentence boundaries set
if "sentencizer" not in nlp_custom.pipe_names:
    print("Adding Sentencizer to custom model...")
    nlp_custom.add_pipe("sentencizer", before="ner")  # Add the sentencizer before the NER component
    print("Sentencizer added successfully.")

# Load terms from a file and return them as a set
def load_terms(filename):
    df = pd.read_csv(filename, sep='\t')
    return set(df['Name'].str.strip().str.lower())

# Load EDC terms and receptor terms
print("Loading EDC terms...")
androgen_terms = load_terms('standardized_EDC_androgen_catalog.tsv')
estrogen_terms = load_terms('standardized_EDC_estrogen_catalog.tsv')
deduct_terms = load_terms('standardized_EDC_catalog_deduct.tsv')
manual_terms = load_terms('manual_catalog.tsv')
edc_terms = androgen_terms.union(estrogen_terms).union(deduct_terms).union(manual_terms)

print("Loading receptor terms...")
receptor_terms = load_terms('standardized_edc_targets.tsv')
additional_receptor_terms = load_terms('receptors.tsv')
receptor_terms.update(additional_receptor_terms)

# Create and save a DataFrame from the union of receptor terms
receptor_terms_df = pd.DataFrame({'Name': list(receptor_terms)})
receptor_terms_df.to_csv('combined_target_catalog.tsv', sep='\t', index=False)
print("The union of receptor terms has been saved to combined_target_catalog.tsv.")

# Function to perform NER and add results to the table
def ner(text, pmid, sentence_id, sentence_map, models):
    model_results = {}  # Track results from each model

    # Process the text with each model
    for model in models:
        doc = model(text.lower())
        abbreviations = {str(abrv).lower(): str(abrv._.long_form).lower() for abrv in doc._.abbreviations}
        for short_form, long_form in abbreviations.items():
            text = text.replace(short_form, long_form)
        doc_expanded = model(text.lower())

        # Record the results from each model
        for ent in doc_expanded.ents:
            ent_text = ent.text.lower()
            if ent_text in edc_terms:
                ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
            elif ent_text in receptor_terms:
                ent_label = "TARGET"
            elif ent.label_ == "ACTIVITY":
                ent_label = "ACTIVITY"
                print(f"ACTIVITY found: {ent_text}")
            else:
                ent_label = "UNKNOWN"  # Default label if not matched

            # Update results based on prioritization
            if model == nlp_edc_custom:
                # Prioritize labels from `nlp_edc_custom` over others
                if ent_text not in model_results or ent_label != "UNKNOWN":
                    model_results[ent_text] = ent_label
            else:
                # For other models, update only if no prior results or if labeled as "ACTIVITY"
                if ent_text not in model_results or ent_label == "ACTIVITY":
                    model_results[ent_text] = ent_label

    # Update sentence_map with the final labels from NER models
    for ent_text, ent_label in model_results.items():
        entity_key = (pmid, sentence_id, ent_text)
        sentence_map[entity_key] = ent_label

# Function to reprocess text for fallback matching
def fallback_labeling(text, pmid, sentence_id, sentence_map):
    doc = nlp_custom(text.lower())  # Use a model to reprocess for fallback
    for ent in doc.ents:
        ent_text = ent.text.lower()
        if ent_text not in sentence_map:
            if ent_text in edc_terms:
                ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
            elif ent_text in receptor_terms:
                ent_label = "TARGET"
            else:
                ent_label = "UNKNOWN"
            # Update sentence_map with fallback labels
            sentence_map[(pmid, sentence_id, ent_text)] = ent_label

# Function to extract entities from the title
def extract_entities_from_title(title, pmid, sentence_map, models):
    text = title.lower()
    ner(text, pmid, "title", sentence_map, models)
    fallback_labeling(text, pmid, "title", sentence_map)  # Apply fallback labeling

# Read the CSV file
print("Reading CSV file...")
try:
    df = pd.read_csv('processed_articles.csv')
except Exception as e:
    print(f"Error reading CSV file: {e}")
    raise

print("CSV file read successfully.")
df = df.fillna("")
df_subset = df.head(200)

# Initialize a dictionary to store unique entities
sentence_map = {}

# Initialize counters
sentences_processed = 0
total_sentences = 0

# Define models list
models = [nlp_bi, nlp_custom, nlp_edc_custom]

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    pmid = row['PMID']
    title = row['Title']
    sentence = row['Sentence']
    sentence_id = row['Sentence id']

    title = str(title)
    sentence = str(sentence)

    print(f"Extracting entities from title for row {index + 1}/{len(df)} - PMID: {pmid}")
    extract_entities_from_title(title, pmid, sentence_map, models)

    doc = nlp_bi(sentence)
    sentences = [sent.text for sent in doc.sents]
    total_sentences += len(sentences)

    for sent_id, sent in enumerate(sentences, start=sentence_id):
        print(f"Processing sentence {sent_id} in abstract for row {index + 1}/{len(df)} - PMID: {pmid}")
        ner(sent, pmid, sent_id, sentence_map, models)
        fallback_labeling(sent, pmid, sent_id, sentence_map)  # Apply fallback labeling
        sentences_processed += 1

print(f"Total sentences processed: {sentences_processed}")
print(f"Total sentences in dataset: {total_sentences}")

# Convert the results dictionary to a DataFrame
entity_df = pd.DataFrame([(pmid, text, label, sent_id) for (pmid, sent_id, text), label in sentence_map.items()],
                         columns=["PMID", "Entity", "Class", "Sentence id"])

entity_df = entity_df.sort_values(by=["PMID", "Sentence id"])
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
