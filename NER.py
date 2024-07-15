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

# Load EDC terms from the combined catalog
print("Loading EDC terms from combined_edc_catalog.tsv...")
edc_df = pd.read_csv('combined_edc_catalog.tsv', sep='\t')
edc_terms = set(edc_df['Name'].str.lower())  # Convert to a set for quick lookup


# Function to perform additional checks for chemical entities
def is_chemical(text):
    exclude_words = ["protein"]
    # Example criteria for identifying chemical entities
    if text.isalpha() and len(text) >= 2 and text.lower() not in exclude_words:
        return True
    return False

# Function to perform NER, standardize abbreviations, and add results to the table
def ner(text, pmcid, table, model_name):
    if model_name == "bi":
        doc = nlp_bi(text)
    elif model_name == "bc":
        doc = nlp_bc(text)
    else:
        print("ERROR: Invalid 'model_name' value")
        return table

    # Process abbreviations
    abbreviations = {}
    for abrv in doc._.abbreviations:
        long_form = str(abrv._.long_form)
        short_form = str(abrv)
        abbreviations[short_form.lower()] = long_form
        text = text.replace(short_form, long_form)
    
    # Lowercase the text
    text = text.lower()
    
    # Reprocess the document with expanded abbreviations
    doc = nlp_bi(text) if model_name == "bi" else nlp_bc(text)

    for ent in doc.ents:
        # Check if the entity is followed by a verb or infinitive verb
        is_followed_by_verb = False
        if ent.end < len(doc):
            next_token = doc[ent.end]
            if next_token.pos_ in ["VERB", "ADV"]: #THE INF NEEDS TO BE ADVERB
                is_followed_by_verb = True
        
        # Check for EDC terms
        if ent.text.lower() in edc_terms:
            ent_label = "ENDOCRINE_DISRUPTING_CHEMICAL"
        else:
            ent_label = ent.label_

 #Filter entities by specific classes
        if ent_label in ["ENDOCRINE_DISRUPTING_CHEMICAL", "SIMPLE_CHEMICAL", "CHEMICAL"]:
            # Add entity to the table
            table["ID"].append(pmcid)
            table["Entity"].append(ent.text)
            table["Class"].append(ent_label)
            table["LongForm"].append(abbreviations.get(ent.text.lower(), ""))
            table["FollowedByVerb"].append(is_followed_by_verb)

    return table

# Read the CSV file
print("Reading CSV file...")
df = pd.read_csv('articles.csv')
print("CSV file read successfully.")

# Replace NaN values with empty strings
df = df.fillna("")

df_subset = df.head(10)  # Adjust the subset size as needed

# Initialize the table for storing entities
table = {"ID": [], "Entity": [], "Class": [], "LongForm": [], "FollowedByVerb": []}

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    pmcid = row['PMC ID']
    title = row['Title']
    abstract = row['Abstract']

    # Ensure the text fields are strings
    title = str(title)
    abstract = str(abstract)

    # Print progress
    print(f"Processing row {index + 1}/{len(df)} - PMC ID: {pmcid}")

    # Apply NER to the title
    ner(title, pmcid, table, "bi")
    ner(title, pmcid, table, "bc")

    # Apply NER to the abstract
    ner(abstract, pmcid, table, "bi")
    ner(abstract, pmcid, table, "bc")

# Convert the results table to a DataFrame and save to a CSV file
entity_df = pd.DataFrame(table)
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
