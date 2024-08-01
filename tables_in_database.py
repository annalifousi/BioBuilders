import pandas as pd

# Load the data from the CSV file
df = pd.read_csv("ner_data.csv")

# Group by 'PM ID' and 'Sentence id', then aggregate the 'Entity' values based on 'Class'
grouped_df = df.groupby(['PM ID', 'Sentence id', 'Class'])['Entity'].apply(lambda x: ', '.join(x)).reset_index()

# Pivot the DataFrame
pivot_df = grouped_df.pivot_table(index=['PM ID', 'Sentence id'], columns='Class', values='Entity', aggfunc='first').reset_index()

# Rename columns for clarity
pivot_df.columns.name = None
pivot_df.columns = ['PM ID', 'Sentence id', 'ACTIVITY', 'ENDOCRINE_DISRUPTING_CHEMICAL', 'TARGET']

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

# Ensure PM IDs are treated as strings
final_df['PM ID'] = final_df['PM ID'].astype(str)

# Group by the relevant columns to count occurrences and collect PMIDs
interaction_summary = final_df.groupby(['ACTIVITY', 'ENDOCRINE_DISRUPTING_CHEMICAL', 'TARGET']).agg(
    counts=('PM ID', 'size'),
    articles=('PM ID', lambda x: ', '.join(x.unique()))
).reset_index()

# Save the DataFrame to a CSV file
interaction_summary.to_csv('final_dataframe.csv', index=False)

print("DataFrame saved to 'final_dataframe.csv'")
