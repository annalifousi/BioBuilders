import pandas as pd

# Load the CSV and TSV files
bioactivity_df = pd.read_csv('edc_bioactivity_data.csv')
edc_catalog_df = pd.read_csv('combined_edc_catalog.tsv', sep='\t')

# Remove rows where 'bioactivity_value' is NaN or equals 0.0
bioactivity_with_value = bioactivity_df.dropna(subset=['bioactivity_value'])
bioactivity_with_value = bioactivity_with_value[bioactivity_with_value['bioactivity_value'] != 0.0]

# Find unique compound names in both datasets
bioactivity_compounds = set(bioactivity_with_value['compound_name'].unique())
edc_catalog_compounds = set(edc_catalog_df['Name'].unique())

# Find the intersection of compounds present in both datasets
common_compounds = bioactivity_compounds.intersection(edc_catalog_compounds)

# Filter the bioactivity data to include only records with compound names present in the combined EDC catalog
common_bioactivity_records = bioactivity_with_value[bioactivity_with_value['compound_name'].isin(common_compounds)]

# Save the filtered data to a new CSV file
common_bioactivity_records.to_csv('filtered_edc_bioactivity_data.csv', index=False)

print("Filtered data saved to 'filtered_edc_bioactivity_data.csv'.")
