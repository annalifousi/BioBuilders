import pandas as pd
import re

# Read the data from the TSV file
df = pd.read_csv('edc_targets.tsv', sep='\t', header=None)

# Rename the first column to 'Name'
df.rename(columns={0: 'Name'}, inplace=True)

# Function to standardize names
def standardize_name(name):
    # Convert to lowercase
    name = name.lower()
    # Remove special characters and extra spaces
    name = re.sub(r'[^a-z0-9\s]', '', name)
    name = re.sub(r'\s+', ' ', name).strip()
    return name

# Apply the function to the relevant column
df['standardized_names'] = df['Name'].apply(standardize_name)

# Select only the standardized names column
standardized_names_df = df[['standardized_names']]

# Write the DataFrame with only standardized names to a new TSV file
standardized_names_df.to_csv('standardized_edc_targets.tsv', sep='\t', index=False, header=False)

print("The standardized names have been saved to 'standardized_edc_targets.tsv'.")
