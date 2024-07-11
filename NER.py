import pandas as pd
import spacy
from scispacy.linking import EntityLinker
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md
from scispacy.abbreviation import AbbreviationDetector

# Load the models
print("Loading models...")
nlp_bi = en_ner_bionlp13cg_md.load()
nlp_bc = en_ner_bc5cdr_md.load()
print("Models loaded successfully.")

# Initialize the Entity linker
print("Initializing EntityLinker...")
# linker = EntityLinker(name="umls", k=30)

# Add EntityLinker and AbbreviationDetector to each model's pipeline
nlp_bi.add_pipe("abbreviation_detector")
#nlp_bi.add_pipe("scispacy_linker", config={"resolve_abbreviations": True, "linker_name": "umls"})
nlp_bc.add_pipe("abbreviation_detector")
#nlp_bc.add_pipe("scispacy_linker", config={"resolve_abbreviations": True, "linker_name": "umls"})

# Add the entity ruler to the pipeline
ruler_bi = nlp_bi.add_pipe("entity_ruler", before="ner")
ruler_bc = nlp_bc.add_pipe("entity_ruler", before="ner")

import pandas as pd

# Read EDCs from androgen_EKDB.tsv file
print("Reading EDCs from androgen_EKDB.tsv...")
androgen_df = pd.read_csv('EDC_androgen_catalog.tsv', sep='\t')
edc_androgen_names = set(androgen_df['Name'].str.lower().unique())  # Get unique EDC names and convert to lowercase

# Read EDCs from estrogen.tsv file
print("Reading EDCs from estrogen.tsv...")
estrogen_df = pd.read_csv('EDC_estrogen_catalog.tsv', sep='\t')
edc_estrogen_names = set(estrogen_df['Name'].str.lower().unique())  # Get unique EDC names and convert to lowercase

# Read EDCs from EDC_catalog_deduct.tsv file
print("Reading EDCs from EDC_catalog_deduct.tsv...")
deduct_df = pd.read_csv('EDC_catalog_deduct.tsv', sep='\t')
edc_deduct_names = set(deduct_df['Name'].str.lower().unique())  # Get unique EDC names and convert to lowercase

print("Reading EDCs from manual_catalog.tsv...")
manual_df = pd.read_csv('manual_catalog.tsv', sep='\t')
manual_df_names = set(manual_df['Name'].str.lower())  # Get unique EDC names and convert to lowercase

# Combine all unique EDC names
all_edc_names = edc_androgen_names.union(edc_estrogen_names).union(edc_deduct_names).union(manual_df_names)


# Convert the set to a sorted list
all_edc_names_sorted = sorted(all_edc_names)

# Save the combined unique entities to a new TSV file
edc_df = pd.DataFrame(all_edc_names_sorted, columns=['Name'])
edc_df.to_csv('combined_edc_catalog.tsv', sep='\t', index=False)

print("Combined EDC catalog created successfully with unique entities.")

patterns = [{"label": "ENDOCRINE_DISRUPTING_CHEMICALS", "pattern": edc} for edc in all_edc_names]

# Add patterns to the entity ruler
ruler_bi.add_patterns(patterns)
ruler_bc.add_patterns(patterns)

from spacy.matcher import PhraseMatcher
from spacy.language import Language
# Initialize the PhraseMatcher


# Initialize the PhraseMatcher
matcher = PhraseMatcher(nlp_bi.vocab)
matcher.add("EDC_TERMS", [nlp_bi.make_doc(term) for term in set(edc_df['Name'].str.lower())])
@Language.component("edc_matcher")
def edc_matcher(doc):
    matches = matcher(doc)
    spans = []
    seen_tokens = set()

    for match_id, start, end in matches:
        span = doc[start:end]
        # Check if any token in the span has already been included in another entity
        if any(token.idx in seen_tokens for token in span):
            continue
        spans.append((start, end))
        seen_tokens.update(token.idx for token in span)

    # Add spans to doc.ents, ensuring no overlapping entities are added
    with doc.retokenize() as retokenizer:
        for start, end in spans:
            retokenizer.merge(doc[start:end])

    return doc # Add the custom component to the pipeline
nlp_bi.add_pipe("edc_matcher", before='ner')
nlp_bc.add_pipe("edc_matcher", before = 'ner')# Function to perform NER, standardize abbreviations, and filter entities

def ner(text, pmcid, table, f):
    try:
        if f == "bi":
            doc = nlp_bi(text)
        elif f == "bc":
            doc = nlp_bc(text)
        else:
            print("ERROR: Invalid 'f' value")
            return table

        # Process abbreviations
        abbreviations = {}
        for abrv in doc._.abbreviations:
            long_form = str(abrv._.long_form)
            short_form = str(abrv)
            text = text.replace(short_form, long_form)
            abbreviations[short_form.lower()] = long_form
        
        # Lowercase the text
        text = text.lower()
        
        # Reprocess the document with expanded abbreviations
        doc = nlp_bi(text) if f == "bi" else nlp_bc(text)

        for ent in doc.ents:
            label = ent.label_.lower()
            
            # Check if the entity has any of the required labels
            if label in ["chemical", "endocrine_disrupting_chemicals", "gene_or_gene_product"]:
                # Determine the primary label based on priority
                if "endocrine_disrupting_chemicals" in label:
                    primary_label = "ENDOCRINE_DISRUPTING_CHEMICALS"
                elif "gene_or_gene_product" in label:
                    primary_label = "GENE_OR_GENE_PRODUCT"
                else:
                    primary_label = ent.label_
                
                # Check if the entity is followed by a noun
                is_followed_by_noun = False
                if ent.end < len(doc):
                    next_token = doc[ent.end]
                    if next_token.pos_ == "NOUN":
                        is_followed_by_noun = True

                if is_chemical(ent.text) and is_followed_by_noun:
                    table["ID"].append(pmcid)
                    table["Entity"].append(ent.text.lower())
                    table["Class"].append(primary_label)
                    # Add long form if abbreviation exists, otherwise empty string
                    if ent.text.lower() in abbreviations:
                        table["LongForm"].append(abbreviations[ent.text.lower()])
                    else:
                        table["LongForm"].append("")
        
        return table
    
    except Exception as e:
        print(f"Error in NER function for pmcid={pmcid}: {str(e)}")
        return table


# Function to perform additional checks for chemical entities
def is_chemical(text):
    exclude_words = ["protein"]
    # Example criteria for identifying chemical entities
    # Adjust based on your specific dataset characteristics
    if text.isalpha() and len(text) >= 2 not in exclude_words:  # Example: must be alphabetic and longer than 2 characters
        return True
    return False
# Read the CSV file
print("Reading CSV file...")
df = pd.read_csv('articles.csv')
print("CSV file read successfully.")

# Replace NaN values with empty strings
df = df.fillna("")

# Process each text in the subset
table = {"ID": [], "Entity": [], "Class": [], "LongForm":[]}
for index, row in df.head(100).iterrows():  # Adjust the subset size as needed
    text = row['Abstract']
    pmcid = row['PMC ID']
    table = ner(text, pmcid, table, "bi")
    table = ner(text, pmcid, table, "bc")

# Convert the table to a DataFrame and save as a TSV file
output_df = pd.DataFrame(table)
output_df.to_csv('filtered_entities.tsv', sep='\t', index=False)