import pandas as pd
import spacy
from scispacy.abbreviation import AbbreviationDetector
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md

# Load the models
print("Loading models...")
nlp_bi = en_ner_bionlp13cg_md.load()
nlp_bc = en_ner_bc5cdr_md.load()
print("Models loaded successfully.")

# Add AbbreviationDetector to each model's pipeline
print("Adding AbbreviationDetector...")
nlp_bi.add_pipe("abbreviation_detector")
nlp_bc.add_pipe("abbreviation_detector")
print("AbbreviationDetector added successfully.")

# Load EDC and receptor terms
def load_terms(filename):
    df = pd.read_csv(filename, sep='\t')
    return set(df['Name'].str.lower())

print("Loading EDC terms from combined_edc_catalog.tsv...")
edc_terms = load_terms('combined_edc_catalog.tsv')

print("Loading receptor terms from receptors.tsv...")
receptor_terms = load_terms('receptors.tsv')

# Function to perform NER and add results to the table
def ner(text, pmid, sentence_id, sentence_map, model):
    # Process document with the given model
    doc = model(text)
    
    # Process abbreviations and create a mapping
    abbreviations = {str(abrv).lower(): str(abrv._.long_form).lower() for abrv in doc._.abbreviations}

    # Expand abbreviations in the text
    expanded_text = text.lower()
    for short_form, long_form in abbreviations.items():
        expanded_text = expanded_text.replace(short_form, long_form)

    # Reprocess the document with expanded abbreviations
    doc_expanded = model(expanded_text)
    
    # Collect entities with sentence IDs
    for sent in doc_expanded.sents:
        for ent in sent.ents:
            # Check for EDC terms and receptor terms
            ent_text = ent.text.lower()
            if ent_text in edc_terms:
                ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
            elif ent_text in receptor_terms:
                ent_label = "RECEPTOR"
            else:
                ent_label = ent.label_

            # Filter entities by specific classes
            if ent_label in ["ENDOCRINE_DISRUPTING_CHEMICAL", "RECEPTOR"]:
                # Add entity to the sentence_map dictionary
                entity_key = (pmid, sentence_id, ent_text)
                if entity_key not in sentence_map:
                    sentence_map[entity_key] = ent_label

# Function to extract entities from the title
def extract_entities_from_title(title, pmid, sentence_map, model):
    # Process the title with the given model
    doc = model(title)
    
    # Process abbreviations and create a mapping
    abbreviations = {str(abrv).lower(): str(abrv._.long_form).lower() for abrv in doc._.abbreviations}

    # Expand abbreviations in the title
    expanded_title = title.lower()
    for short_form, long_form in abbreviations.items():
        expanded_title = expanded_title.replace(short_form, long_form)

    # Reprocess the title with expanded abbreviations
    doc_expanded = model(expanded_title)
    
    # Check for entities in the title
    for ent in doc_expanded.ents:
        ent_text = ent.text.lower()
        if ent_text in edc_terms:
            ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
        elif ent_text in receptor_terms:
            ent_label = "RECEPTOR"
        else:
            continue

        # Add only the first occurrence of each entity to the sentence_map
        entity_key = (pmid, "title", ent_text)
        if entity_key not in sentence_map:
            sentence_map[entity_key] = ent_label

# Read the CSV file
print("Reading CSV file...")
try:
    df = pd.read_csv('processed_articles.csv')
except Exception as e:
    print(f"Error reading CSV file: {e}")
    raise

print("CSV file read successfully.")

# Replace NaN values with empty strings
df = df.fillna("")

# Subset the dataframe for testing (adjust the subset size as needed)
df_subset = df.head(500)

# Initialize a dictionary to store unique entities
sentence_map = {}

# Initialize counters
sentences_processed = 0
total_sentences = 0

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    pmid = row['PMID']
    title = row['Title']
    sentence = row['Sentence']
    sentence_id = row['Sentence id']

    # Ensure the text fields are strings
    title = str(title)
    sentence = str(sentence)

    # Extract entities from the title
    print(f"Extracting entities from title for row {index + 1}/{len(df_subset)} - PM ID: {pmid}")
    extract_entities_from_title(title, pmid, sentence_map, nlp_bi)
    extract_entities_from_title(title, pmid, sentence_map, nlp_bc)
    
    # Process the abstract
    doc = nlp_bi(sentence)  # Use SpaCy's built-in sentence segmentation
    sentences = [sent.text for sent in doc.sents]

    # Track total sentences
    total_sentences += len(sentences)
    
    # Process each sentence individually
    for sent_id, sent in enumerate(sentences, start=sentence_id):
        # Print progress
        print(f"Processing sentence {sent_id} in abstract for row {index + 1}/{len(df_subset)} - PM ID: {pmid}")

        # Apply NER to each sentence using both models
        ner(sent, pmid, sent_id, sentence_map, nlp_bi)
        ner(sent, pmid, sent_id, sentence_map, nlp_bc)
        
        # Increment the counter for each sentence processed
        sentences_processed += 1

print(f"Total sentences processed: {sentences_processed}")
print(f"Total sentences in dataset: {total_sentences}")

# Convert the results dictionary to a DataFrame
entity_df = pd.DataFrame([(pmid, text, label, sent_id) for (pmid, sent_id, text), label in sentence_map.items()],
                         columns=["PM ID", "Entity", "Class", "Sentence id"])

# Sort the DataFrame by PM ID and Sentence id
entity_df = entity_df.sort_values(by=["PM ID", "Sentence id"])

# Save the sorted DataFrame to a CSV file
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
