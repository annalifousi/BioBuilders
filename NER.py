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
def ner(text, pmid, table, model):
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
    for sent_id, sentence in enumerate(doc_expanded.sents):
        for ent in sentence.ents:
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
                # Add entity to the table with sentence ID
                table.add((pmid, ent.text, ent_label, sent_id))

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

# Initialize a set for storing unique entities
table = set()

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    pmid = row['PMID']
    title = row['Title']
    abstract = row['Sentence']
    
    # Ensure the text fields are strings
    title = str(title)
    abstract = str(abstract)

    # Combine title and abstract for processing
    combined_text = f"{title}. {abstract}"

    # Print progress
    print(f"Processing row {index + 1}/{len(df)} - PM ID: {pmid}")

    # Apply NER to the combined text using both models
    ner(combined_text, pmid, table, nlp_bi)
    ner(combined_text, pmid, table, nlp_bc)

# Convert the results set to a DataFrame
entity_df = pd.DataFrame(table, columns=["PM ID", "Entity", "Class", "Sentence id"])

# Sort the DataFrame by PM ID and Sentence id
entity_df = entity_df.sort_values(by=["PM ID", "Sentence id"])

# Save the sorted DataFrame to a CSV file
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
