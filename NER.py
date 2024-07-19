import pandas as pd
import spacy
from scispacy.abbreviation import AbbreviationDetector
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md
from sklearn.feature_extraction.text import CountVectorizer

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

# Load EDC terms from the combined catalog
print("Loading EDC terms from combined_edc_catalog.tsv...")
edc_df = pd.read_csv('combined_edc_catalog.tsv', sep='\t')
edc_terms = set(edc_df['Name'].str.lower())  # Convert to a set for quick lookup

# Load receptor terms from the receptor catalog
print("Loading receptor terms from receptors.tsv...")
receptor_df = pd.read_csv('receptors.tsv', sep='\t')
receptor_terms = set(receptor_df['Name'].str.lower())  # Convert to a set for quick lookup

# Function to perform NER, standardize abbreviations, and add results to the table
def ner(text, pmcid, table, model):
    # Process document with the given model
    doc = model(text)

    # Process abbreviations
    abbreviations = {str(abrv).lower(): str(abrv._.long_form) for abrv in doc._.abbreviations}
    for short_form, long_form in abbreviations.items():
        text = text.replace(short_form, long_form)

    # Lowercase the text
    text = text.lower()

    # Reprocess the document with expanded abbreviations
    doc = model(text)

    # Collect entities with sentence IDs
    for sent_id, sentence in enumerate(doc.sents):
        for ent in sentence.ents:
            # Check for EDC terms
            if ent.text.lower() in edc_terms:
                ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
            elif ent.text.lower() in receptor_terms:
                ent_label = "RECEPTOR"
            else:
                ent_label = ent.label_

            # Filter entities by specific classes
            if ent_label in ["ENDOCRINE_DISRUPTING_CHEMICAL", "RECEPTOR"]:
                # Add entity to the table with sentence ID
                table.add((pmcid, ent.text, ent_label, sent_id))

# Read the CSV file
print("Reading CSV file...")
df = pd.read_csv('processed_articles.csv')
print("CSV file read successfully.")

# Replace NaN values with empty strings
df = df.fillna("")

# Subset the dataframe for testing (adjust the subset size as needed)
df_subset = df.head(500)

# Initialize a set for storing unique entities
table = set()

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    pmcid = row['PMCID']
    title = row['Title']
    sentence = row['Sentence']
    
    # Ensure the text fields are strings
    title = str(title)
    sentence = str(sentence)

    # Print progress
    print(f"Processing row {index + 1}/{len(df)} - PMC ID: {pmcid}")

    # Apply NER to the title and sentence using both models
    for text in [title, sentence]:
        ner(text, pmcid, table, nlp_bi)
        ner(text, pmcid, table, nlp_bc)

# Convert the results set to a DataFrame and save to a CSV file
entity_df = pd.DataFrame(table, columns=["PMC ID", "Entity", "Class", "Sentence id"])
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
