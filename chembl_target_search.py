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
activity = new_client.activity


# Function to search for bioactivity and targets by compound name, filtering by Homo sapiens
def search_bioactivity_by_compound(compound_name):
    print(f"Searching for bioactivity data of compound: {compound_name}")
    compound_data = []
    try:
        mols = molecule.filter(pref_name__iexact=compound_name)
        if mols:
            mol_id = mols[0]['molecule_chembl_id']

            # Retrieve bioactivity data
            activities = activity.filter(molecule_chembl_id=mol_id).only(
                ['target_chembl_id', 'standard_type', 'standard_value', 'standard_units'])

            for act in activities:
                target_id = act['target_chembl_id']
                target_info = target.get(target_id)

                # Filter for targets in Homo sapiens
                if target_info and target_info.get('organism') == 'Homo sapiens':
                    bioactivity_info = {
                        'compound_name': compound_name,
                        'target_name': target_info.get('pref_name', 'Unknown'),
                        'bioactivity_type': act['standard_type'],
                        'bioactivity_value': act['standard_value'],
                        'bioactivity_units': act['standard_units']
                    }
                    compound_data.append(bioactivity_info)

            print(
                f"Found {len(compound_data)} bioactivity records for Homo sapiens targets of compound: {compound_name}")
            return compound_data
        else:
            print(f"No molecule found for compound: {compound_name}")
            return []
    except Exception as e:
        print(f"Error searching for compound {compound_name}: {e}")
        return []


# Iterate over the compounds and search for bioactivity data
all_bioactivity_data = []
for i, compound_name in enumerate(compound_names):
    print(f"Processing {i + 1}/{len(compound_names)}: {compound_name}")
    compound_bioactivity = search_bioactivity_by_compound(compound_name)
    all_bioactivity_data.extend(compound_bioactivity)

# Convert to DataFrame
bioactivity_df = pd.DataFrame(all_bioactivity_data)

# Save bioactivity data to CSV file
output_file_path = 'edc_bioactivity_data.csv'  # Replace with your desired output file path
bioactivity_df.to_csv(output_file_path, index=False)
print(f"Bioactivity data saved to {output_file_path}")
