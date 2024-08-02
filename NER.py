import pandas as pd
import spacy
from spacy.pipeline import Sentencizer
from scispacy.abbreviation import AbbreviationDetector
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md
from transformers import BertTokenizer, BertForTokenClassification, pipeline

# Load SpaCy models
print("Loading models...")
nlp_bi = en_ner_bionlp13cg_md.load()
nlp_bc = en_ner_bc5cdr_md.load()
nlp_custom = spacy.load('./model-best')
print("Models loaded successfully.")

# Add abbreviation detector
print("Adding AbbreviationDetector...")
nlp_bi.add_pipe("abbreviation_detector")
nlp_bc.add_pipe("abbreviation_detector")
nlp_custom.add_pipe("abbreviation_detector")
print("AbbreviationDetector added successfully.")

# Add sentencizer to custom model
if "sentencizer" not in nlp_custom.pipe_names:
    print("Adding Sentencizer to custom model...")
    nlp_custom.add_pipe("sentencizer", before="ner")
    print("Sentencizer added successfully.")

# Load EDC and receptor terms
def load_terms(filename):
    df = pd.read_csv(filename, sep='\t')
    return set(df['Name'].str.lower())

print("Loading EDC and receptor terms...")
edc_terms = load_terms('combined_edc_catalog.tsv')
receptor_terms = load_terms('standardized_edc_targets.tsv')
additional_receptor_terms = load_terms('receptors.tsv')
receptor_terms.update(additional_receptor_terms)

# Load SciBERT model for contextual embeddings
print("Loading SciBERT model...")
tokenizer = BertTokenizer.from_pretrained('allenai/scibert_scivocab_cased')
model = BertForTokenClassification.from_pretrained('allenai/scibert_scivocab_cased')

# Load pipeline for NER using SciBERT
nlp_sci_bert = pipeline('ner', model=model, tokenizer=tokenizer)

# Initialize SpaCy lemmatizer
print("Initializing SpaCy lemmatizer...")
lemmatizer = spacy.load('en_core_web_sm')

# Function for contextual entity expansion
def contextual_expansion(text, model):
    expanded_entities = set()
    results = model(text)
    for result in results:
        if result['entity'] in ['LABEL_1', 'LABEL_2']:  # Adjust based on the entity labels from SciBERT
            expanded_entities.add(result['word'].lower())
    return expanded_entities

# Function to lemmatize activity entities
def lemmatize_text(text):
    doc = lemmatizer(text)
    return " ".join([token.lemma_ for token in doc])

# Function to perform NER and add results to the table
def ner(text, pmid, sentence_id, sentence_map, model, contextual_model):
    doc = model(text)
    abbreviations = {str(abrv).lower(): str(abrv._.long_form).lower() for abrv in doc._.abbreviations}
    expanded_text = text.lower()
    for short_form, long_form in abbreviations.items():
        expanded_text = expanded_text.replace(short_form, long_form)
    
    # Contextual expansion
    expanded_entities = contextual_expansion(expanded_text, contextual_model)
    
    doc_expanded = model(expanded_text)
    for sent in doc_expanded.sents:
        for ent in sent.ents:
            ent_text = ent.text.lower()
            is_lemmatized_activity = False
            if ent_text in edc_terms or ent_text in expanded_entities:
                ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
            elif ent_text in receptor_terms:
                ent_label = "TARGET"
            elif ent.label_ == "ACTIVITY":
                ent_label = "ACTIVITY"
                # Lemmatize the activity term
                lemmatized_text = lemmatize_text(ent.text)
                print(f"ACTIVITY found: {ent_text} (lemmatized: {lemmatized_text})")
                ent_text = lemmatized_text  # Use lemmatized text for consistency
                is_lemmatized_activity = True
            else:
                continue

            if ent_label in ["ENDOCRINE_DISRUPTING_CHEMICAL", "TARGET", "ACTIVITY"]:
                entity_key = (pmid, sentence_id, ent_text)
                if entity_key not in sentence_map:
                    sentence_map[entity_key] = {
                        "label": ent_label,
                        "is_lemmatized_activity": is_lemmatized_activity
                    }

# Function to extract entities from the title
def extract_entities_from_title(title, pmid, sentence_map, model, contextual_model):
    doc = model(title)
    abbreviations = {str(abrv).lower(): str(abrv._.long_form).lower() for abrv in doc._.abbreviations}
    expanded_title = title.lower()
    for short_form, long_form in abbreviations.items():
        expanded_title = expanded_title.replace(short_form, long_form)

    # Contextual expansion
    expanded_entities = contextual_expansion(expanded_title, contextual_model)
    
    doc_expanded = model(expanded_title)
    for ent in doc_expanded.ents:
        ent_text = ent.text.lower()
        is_lemmatized_activity = False
        if ent_text in edc_terms or ent_text in expanded_entities:
            ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
        elif ent_text in receptor_terms:
            ent_label = "TARGET"
        elif ent.label_ == "ACTIVITY":
            ent_label = "ACTIVITY"
            # Lemmatize the activity term
            lemmatized_text = lemmatize_text(ent.text)
            print(f"ACTIVITY found in title: {ent_text} (lemmatized: {lemmatized_text})")
            ent_text = lemmatized_text  # Use lemmatized text for consistency
            is_lemmatized_activity = True
        else:
            continue

        entity_key = (pmid, "title", ent_text)
        if entity_key not in sentence_map:
            sentence_map[entity_key] = {
                "label": ent_label,
                "is_lemmatized_activity": is_lemmatized_activity
            }

# Main processing
print("Reading CSV file...")
try:
    df = pd.read_csv('processed_articles.csv')
except Exception as e:
    print(f"Error reading CSV file: {e}")
    raise

print("CSV file read successfully.")
df = df.fillna("")
df_subset = df.head(500)
sentence_map = {}
sentences_processed = 0
total_sentences = 0

for index, row in df_subset.iterrows():
    pmid = row['PMID']
    title = str(row['Title'])
    sentence = str(row['Sentence'])
    sentence_id = row['Sentence id']

    print(f"Extracting entities from title for row {index + 1}/{len(df)} - PM ID: {pmid}")
    extract_entities_from_title(title, pmid, sentence_map, nlp_bi, nlp_sci_bert)
    extract_entities_from_title(title, pmid, sentence_map, nlp_bc, nlp_sci_bert)
    extract_entities_from_title(title, pmid, sentence_map, nlp_custom, nlp_sci_bert)

    doc = nlp_bi(sentence)
    sentences = [sent.text for sent in doc.sents]
    total_sentences += len(sentences)

    for sent_id, sent in enumerate(sentences, start=sentence_id):
        print(f"Processing sentence {sent_id} in abstract for row {index + 1}/{len(df)} - PM ID: {pmid}")
        ner(sent, pmid, sent_id, sentence_map, nlp_bi, nlp_sci_bert)
        ner(sent, pmid, sent_id, sentence_map, nlp_bc, nlp_sci_bert)
        ner(sent, pmid, sent_id, sentence_map, nlp_custom, nlp_sci_bert)

        sentences_processed += 1

print(f"Total sentences processed: {sentences_processed}")
print(f"Total sentences in dataset: {total_sentences}")

# Convert sentence_map to DataFrame
entity_df = pd.DataFrame([
    (pmid, text, info["label"], sent_id, info["is_lemmatized_activity"])
    for (pmid, sent_id, text), info in sentence_map.items()
], columns=["PM ID", "Entity", "Class", "Sentence id", "Is Lemmatized Activity"])

entity_df = entity_df.sort_values(by=["PM ID", "Sentence id"])
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
