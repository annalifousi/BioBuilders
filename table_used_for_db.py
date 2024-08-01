import pandas as pd

# Load the data
file_path = 'ner_data.csv'
df = pd.read_csv(file_path, encoding='ascii')

import pandas as pd

# Load the data
file_path = 'ner_data.csv'
df = pd.read_csv(file_path, encoding='ascii')

# Filter for ENDOCRINE_DISRUPTING_CHEMICAL and RECEPTOR classes
edc_df = df[df['Class'] == 'ENDOCRINE_DISRUPTING_CHEMICAL']
receptor_df = df[df['Class'] == 'RECEPTOR']

# Merge on PMC ID and Sentence id to find associations
merged_df = pd.merge(edc_df, receptor_df, on=['PMC ID', 'Sentence id'], suffixes=('_EDC', '_RECEPTOR'))

# Display the merged dataframe
print(merged_df.head())