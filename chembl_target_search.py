import pandas as pd
from chembl_webresource_client.new_client import new_client

# Load the TSV file
tsv_file_path = 'combined_edc_catalog.tsv'  # Replace with your file path
print(f"Loading TSV file from: {tsv_file_path}")
edc_df = pd.read_csv(tsv_file_path, sep='\t')

# Print the column names to inspect them
print("Column names in the TSV file:", edc_df.columns)

# Assuming the actual column name is 'Name' (replace it with the correct name)
compound_column_name = 'Name'  # Adjust this to the correct column name in your file

# Extract compound names
compound_names = edc_df[compound_column_name].tolist()
print(f"Extracted {len(compound_names)} compound names.")

# Initialize ChEMBL clients
molecule = new_client.molecule
target = new_client.target

# Function to search for targets by compound name
def search_targets_by_compound(compound_name):
    print(f"Searching for targets of compound: {compound_name}")
    try:
        mols = molecule.filter(pref_name__iexact=compound_name)
        if mols:
            mol_id = mols[0]['molecule_chembl_id']
            activities = new_client.activity.filter(molecule_chembl_id=mol_id).only(['target_chembl_id'])
            target_ids = set(act['target_chembl_id'] for act in activities)
            targets = [target.get(tid) for tid in target_ids]
            human_targets = [t for t in targets if t and t.get('organism') == 'Homo sapiens']
            print(f"Found {len(human_targets)} human targets for compound: {compound_name}")
            return human_targets
        else:
            print(f"No molecule found for compound: {compound_name}")
            return []
    except Exception as e:
        print(f"Error searching for compound {compound_name}: {e}")
        return []

# Iterate over the first 100 compounds and search for targets
all_targets = {}
for i, compound_name in enumerate(compound_names):
    print(f"Processing {i+1}/{len(compound_names)}: {compound_name}")
    targets = search_targets_by_compound(compound_name)
    all_targets[compound_name] = targets

# Extract only target names for saving
print("Saving results to file...")
output_file_path = 'edc_targets.tsv'  # Replace with your desired output file path
with open(output_file_path, 'w') as f:
    for compound_name, targets in all_targets.items():
        for target in targets:
            if target and 'pref_name' in target:
                f.write(f"{target['pref_name']}\n")
print("Process completed successfully.")


targets = pd.read_csv('edc_targets.tsv', sep='\t')

unique_targets = targets.drop_duplicates()
unique_targets.to_csv('edc_targets.tsv', sep='\t', index=False)
