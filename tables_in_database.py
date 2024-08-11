import pandas as pd

# Load the data from the CSV file
df = pd.read_csv("updated_filtered_entity_rec.csv")

# Filter the DataFrame to keep only the specified classes
df = df[df['Label'].isin(['ENDOCRINE_DISRUPTING_CHEMICAL', 'TARGET', 'ACTIVITY'])]

# Group by 'PMID' and 'SentenceID', then aggregate the 'Entity' values by their label
grouped_df = df.groupby(['PMID', 'SentenceID', 'Label'])['Entity'].apply(lambda x: ', '.join(x)).reset_index()

# Pivot the DataFrame to have each label as a separate column
pivot_df = grouped_df.pivot_table(index=['PMID', 'SentenceID'], columns='Label', values='Entity', aggfunc='first').reset_index()

# Rename columns for clarity
pivot_df.columns.name = None
pivot_df.columns = ['PMID', 'SentenceID', 'ACTIVITY', 'ENDOCRINE_DISRUPTING_CHEMICAL', 'TARGET']

# Fill NaN values with empty strings
pivot_df.fillna('', inplace=True)

# Function to split columns into separate rows
def split_column_to_rows(df, column):
    rows = []
    for _, row in df.iterrows():
        entities = row[column].split(', ')
        for entity in entities:
            if entity:
                new_row = row.copy()
                new_row[column] = entity
                rows.append(new_row)
    return pd.DataFrame(rows)

# Apply the splitting function to each column
activity_df = split_column_to_rows(pivot_df, 'ACTIVITY')
endocrine_df = split_column_to_rows(activity_df, 'ENDOCRINE_DISRUPTING_CHEMICAL')
final_df = split_column_to_rows(endocrine_df, 'TARGET')

# Reset index for clean output
final_df = final_df.reset_index(drop=True)

# Ensure PMIDs are treated as strings
final_df['PMID'] = final_df['PMID'].astype(str)

# Group by the relevant columns to count occurrences and collect PMIDs
interaction_summary = final_df.groupby(['ACTIVITY', 'ENDOCRINE_DISRUPTING_CHEMICAL', 'TARGET']).agg(
    counts=('PMID', 'size'),
    articles=('PMID', lambda x: ', '.join(x.unique()))
).reset_index()

# Reorder the columns
interaction_summary = interaction_summary[['ENDOCRINE_DISRUPTING_CHEMICAL', 'ACTIVITY', 'TARGET', 'counts', 'articles']]

# Save the DataFrame to a CSV file
interaction_summary.to_csv('final_dataframe.csv', index=False)

print("DataFrame saved to 'final_dataframe.csv'")
