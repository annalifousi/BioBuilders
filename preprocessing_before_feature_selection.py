import pandas as pd

# Function to load terms from the catalogs
def load_terms(filename):
    df = pd.read_csv(filename, sep='\t')
    return set(df['Name'].str.lower())

# Load catalogs
edc_terms = load_terms('combined_edc_catalog.tsv')
receptor_terms = load_terms('standardized_edc_targets.tsv')
additional_receptor_terms = load_terms('receptors.tsv')
receptor_terms.update(additional_receptor_terms)

# EDC name alternatives
edc_alternatives = set([
    'endocrine disrupting chemicals', 
    'sex hormone disruptors', 
    'endocrine disrupting compounds', 
    'endocrine-disrupting chemicals', 
    'endocrine-disrupting compounds', 
    'hormone antagonist', 
    'sex hormone antagonist'
])
# Combine EDC terms with EDC alternatives
edc_terms.update(edc_alternatives)

# Load the main data
df = pd.read_csv('final_annotated_entities.csv')

# Convert Entity column to lowercase
df['Entity'] = df['Entity'].str.lower()

# Ensure all entities are strings
df['Entity'] = df['Entity'].astype(str)

# Function to check for matches in the catalogs
def match_entity(entities):
    combined_entity = ' '.join(entities)
    if combined_entity in edc_terms:
        return combined_entity, 'ENDOCRINE_DISRUPTING_CHEMICAL'
    elif combined_entity in receptor_terms:
        return combined_entity, 'TARGET'
    return None, None

# Create a new DataFrame to store updated rows
updated_rows = []

# Iterate over each group of rows with the same PMID and SentenceID
for _, group in df.groupby(['PMID', 'SentenceID']):
    entities = group['Entity'].tolist()
    labels = group['Label'].tolist()
    i = 0
    
    while i < len(entities):
        matched = False
        
        # Try to match in steps of 4, 3, or 2 words
        for step in range(4, 1, -1):
            if i + step <= len(entities):
                candidate_entities = entities[i:i + step]
                combined_entity, label = match_entity(candidate_entities)
                if label:
                    # Combine original labels if they exist (to avoid losing them)
                    for j in range(step):
                        if labels[i + j] != 'UNKNOWN':
                            label = labels[i + j]  # Use the original label if it's not 'UNKNOWN'
                            break
                    
                    updated_rows.append([group['PMID'].iloc[0], group['SentenceID'].iloc[0], combined_entity, label])
                    i += step  # Skip the matched entities
                    matched = True
                    break
        
        # If no match is found, keep the entity as is, but combine with any existing label
        if not matched:
            updated_rows.append([group['PMID'].iloc[0], group['SentenceID'].iloc[0], entities[i], labels[i]])
            i += 1

# Create a new DataFrame with the updated rows
updated_df = pd.DataFrame(updated_rows, columns=['PMID', 'SentenceID', 'Entity', 'Label'])

# Save the updated DataFrame to a new CSV file
updated_df.to_csv('updated_filtered_entity_rec.csv', index=False)