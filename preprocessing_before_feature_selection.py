import pandas as pd
import inflect
import string

#This steps belongs to the post processing after we have conducted the NER 

#account for the plurals also in the tsv files
def get_plurals(word):
    plurals = set()
    plurals.add(word)
    plural = p.plural(word)
    if plural:
        plurals.add(plural)
    return plurals

#load the tsv files 
def load_terms(filename):
    try:
        df = pd.read_csv(filename, sep='\t')
        terms = df['Name'].str.lower().str.strip()
        terms_with_plurals = set()
        for term in terms:
            for plural in get_plurals(term):
                terms_with_plurals.add(plural)
        logging.debug(f"Loaded terms from {filename}: {terms_with_plurals}")
        return terms_with_plurals
    except Exception as e:
        logging.error(f"Error loading terms from {filename}: {e}")
        raise

# Load catalogs
edc_terms = load_terms('combined_edc_catalog.tsv')
receptor_terms = load_terms('combined_target_catalog.tsv')

# Define EDC name alternatives
edc_alternatives = set([
    'endocrine disrupting chemicals', 'sex hormone disruptors', 'endocrine disrupting compounds', 
    'endocrine-disrupting chemicals', 'endocrine-disrupting compounds', 'hormone antagonist', 
    'sex hormone antagonist', 'endocrine-disrupting contaminants',
])

# Define target alternatives
target_alternatives = set([
    'steroid receptors', 'steroid receptor', 'androgen receptors', 'estrogen receptors', 
    'hormone receptors', 'hormone receptor'
])

# Combine terms and alternatives
edc_terms.update(edc_alternatives)
receptor_terms.update(target_alternatives)

# Load the main data
df = pd.read_csv('ner_data.csv')
df['Entity'] = df['Entity'].str.lower().str.strip().astype(str)

def clean_entity(entity):
    cleaned = entity.strip().rstrip(string.punctuation)
    return cleaned

def match_entity(entity):
    cleaned_entity = clean_entity(entity)
    if cleaned_entity in edc_terms:
        return 'ENDOCRINE_DISRUPTING_CHEMICAL'
    if cleaned_entity in receptor_terms:
        return 'TARGET'
    return None

def additional_conditions(entity):
    words = entity.split()
    if not words:
        return None

    last_word = words[-1]
    if last_word in ['receptor', 'receptors']:
        if len(words) > 1 and words[-2] in ['agonist', 'antagonist', 'agonists', 'antagonists', 'ligand', 'ligands']:
            return 'ENDOCRINE_DISRUPTING_CHEMICAL'
        return 'TARGET'
    elif last_word in ['agonist', 'antagonist', 'agonists', 'antagonists', 'ligand', 'ligands','inhibitor','inhibitors','binder','binders]':
        return 'ENDOCRINE_DISRUPTING_CHEMICAL'
    return None

def process_group(group):
    entities = group['Entity'].tolist()
    updated_rows = []

    for window_size in range(4, 0, -1):
        for i in range(len(entities) - window_size + 1):
            combined_entity = ' '.join(entities[i:i + window_size]).strip()

            # Attempt to match the entity
            new_label = match_entity(combined_entity)
            if new_label:
                final_label = new_label
                additional_label = additional_conditions(combined_entity)
                if additional_label:
                    final_label = additional_label
                updated_rows.append([group['PMID'].iloc[0], group['Sentence id'].iloc[0], combined_entity, final_label])
                break
        else:
            # If no match was found, keep the original label
            updated_rows = group[['PMID', 'Sentence id', 'Entity', 'Class']].values.tolist()
    
    return updated_rows

# Process data
activity_rows = df[df['Class'] == 'ACTIVITY']
updated_rows = []

for _, group in df[df['Class'] != 'ACTIVITY'].groupby(['PMID', 'Sentence id']):
    try:
        updated_rows.extend(process_group(group))
    except Exception as e:
        logging.error(f"Error processing group with PMID {group['PMID'].iloc[0]} and SentenceID {group['Sentence id'].iloc[0]}: {e}")

# Create the updated DataFrame
updated_df = pd.DataFrame(updated_rows, columns=['PMID', 'Sentence id', 'Entity', 'Class'])
final_df = pd.concat([updated_df, activity_rows], ignore_index=True)
final_df_sorted = final_df.sort_values(by=['PMID', 'Sentence id']).reset_index(drop=True)
final_df_sorted.to_csv('updated_ner_data.csv', index=False)

