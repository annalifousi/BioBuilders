import pandas as pd
import spacy
from spacy.pipeline import Sentencizer
from scispacy.abbreviation import AbbreviationDetector
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md
import zipfile
import os

import spacy
import zipfile
import os

# Function to unzip and load a spaCy model
def load_model_from_zip(zip_file_path):
    # Unzip the file to a temporary directory
    model_dir = zip_file_path.replace('.zip', '')  # Create a directory name without the '.zip'
    
    # Check if the directory already exists (from a previous extraction)
    if not os.path.exists(model_dir):
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(model_dir)
    
    # Check if meta.json exists in the immediate directory
    meta_json_path = os.path.join(model_dir, 'meta.json')
    
    if not os.path.exists(meta_json_path):
        # If meta.json is not found, check subdirectories
        for root, dirs, files in os.walk(model_dir):
            if 'meta.json' in files:
                model_dir = root
                break
        else:
            raise IOError(f"Could not find meta.json in {zip_file_path}")

    # Load the model using spaCy
    nlp_model = spacy.load(model_dir)
    return nlp_model

print("Loading models...")

# Load the pretrained spaCy models
nlp_bi = spacy.load("en_ner_bionlp13cg_md")  # Pretrained SCISPACY NER model for biomedical entities
nlp_bc = spacy.load("en_ner_bc5cdr_md")      # Pretrained SCISPACY NER model for biomedical entities

# Load the custom models from zip files
nlp_custom = load_model_from_zip('./model-best.zip')          # Custom model for entities that reveal ACTIVITY
nlp_edc_custom = load_model_from_zip('./EDC_model-best.zip')  # Refined nlp_bi model to recognize EDCs and their molecular TARGETS

print("Models loaded successfully.")



# Add AbbreviationDetector to each model's pipeline to capture the abreviations and the long form of the entities
print("Adding AbbreviationDetector...")
for nlp in [nlp_bi, nlp_bc, nlp_custom]:
    nlp.add_pipe("abbreviation_detector")
print("AbbreviationDetector added successfully.")

# Add a sentencizer to the custom model if it does not have sentence boundaries set
if "sentencizer" not in nlp_custom.pipe_names:
    print("Adding Sentencizer to custom model...")
    nlp_custom.add_pipe("sentencizer", before="ner")  # Add the sentencizer before the NER component
    print("Sentencizer added successfully.")

# Function to perform NER and label the terms appropriately
def ner(text, pmid, sentence_id, sentence_map, models):
    # Process the text with each model
    for model in models:
        doc = model(text.lower())
        abbreviations = {str(abrv).lower(): str(abrv._.long_form).lower() for abrv in doc._.abbreviations}
        for short_form, long_form in abbreviations.items():
            text = text.replace(short_form, long_form)
        doc_expanded = model(text.lower())

        # Debugging: Print detected entities
        print(f"Entities detected by {model}: {[ent.text.lower() for ent in doc_expanded.ents]}")

        # Record the results from each model
        for ent in doc_expanded.ents:
            ent_text = ent.text.lower()
            ent_label = ent.label_
            sentence_map[(pmid, sentence_id, ent_text)] = ent_label

    # Annotate all terms with "UNKNOWN" if not already labeled. We will have more words in the ner_data.csv because we keep all the keywords but this will help us with the post processing the filtering and the curation of the results of the NER
    doc_full = spacy.tokens.Doc(nlp_bi.vocab, words=[token.text.lower() for token in nlp_bi(text.lower()).doc])
    for token in doc_full:
        term = token.text.lower()
        if (pmid, sentence_id, term) not in sentence_map:
            sentence_map[(pmid, sentence_id, term)] = "UNKNOWN"

# Function to extract entities from the title
def extract_entities_from_title(title, pmid, sentence_map, models):
    text = title.lower()
    ner(text, pmid, "title", sentence_map, models)

# Read the CSV file
print("Reading CSV file...")
try:
    df = pd.read_csv('processed_articles.csv')
except Exception as e:
    print(f"Error reading CSV file: {e}")
    raise

print("CSV file read successfully.")
df = df.fillna("")
df_subset = df.head(200) # this is the trial dataset.

# Initialize a dictionary to store unique entities
sentence_map = {}

# Define models list that we are going to perform the NER
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

    for sent_id, sent in enumerate(sentences, start=sentence_id):
        print(f"Processing sentence {sent_id} in abstract for row {index + 1}/{len(df)} - PMID: {pmid}")
        ner(sent, pmid, sent_id, sentence_map, models)

print(f"Total sentences processed: {len(sentence_map)}")

# Convert the results dictionary to a DataFrame
entity_df = pd.DataFrame([(pmid, text, label, sent_id) for (pmid, sent_id, text), label in sentence_map.items()],
                         columns=["PMID", "Entity", "Class", "Sentence id"])

entity_df = entity_df.sort_values(by=["PMID", "Sentence id"])
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
