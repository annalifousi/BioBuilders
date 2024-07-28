import pandas as pd

# Load the CSV files
csv1 = pd.read_csv('ner_data.csv')
csv2 = pd.read_csv('processed_articles.csv')

# Ensure the Sentence id and PM ID columns are of the same type (convert to string)
csv1['Sentence id'] = csv1['Sentence id'].astype(str)
csv2['Sentence id'] = csv2['Sentence id'].astype(str)
csv1['PM ID'] = csv1['PM ID'].astype(str)
csv2['PMID'] = csv2['PMID'].astype(str)

# Define the classes of interest
classes_of_interest = ['ENDOCRINE_DISRUPTING_CHEMICAL', 'ACTIVITY', 'TARGET']

# Filter csv1 for entities in the classes of interest
filtered_csv1 = csv1[csv1['Class'].isin(classes_of_interest)]

# Group by PM ID and Sentence id and filter groups that contain all three classes
grouped = filtered_csv1.groupby(['PM ID', 'Sentence id']).filter(
    lambda x: set(classes_of_interest).issubset(set(x['Class']))
)

# Merge the filtered groups with the second CSV on PM ID and Sentence id
merged_df = pd.merge(grouped, csv2, left_on=['PM ID', 'Sentence id'], right_on=['PMID', 'Sentence id'])

# Select the required columns
result_df = merged_df[['PMID', 'Sentence id', 'Sentence', 'Entity', 'Class']]

# Count unique sentences that meet the criteria
unique_sentences_count = result_df[['PMID', 'Sentence id']].drop_duplicates().shape[0]

# Save the result to a new CSV file
output_path = 'filtered_sentences.csv'  # Update with your desired output path
result_df.to_csv(output_path, index=False)

print("Filtered results saved to:", output_path)
print(f"Number of unique sentences that meet the criteria: {unique_sentences_count}")