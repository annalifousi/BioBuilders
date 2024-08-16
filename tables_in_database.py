import pandas as pd
from itertools import product

# Load the data from the CSV file
df = pd.read_csv('filtered_sentences.csv', sep=',')

# Pivot table to get the desired format
pivot_table = pd.pivot_table(df, 
                             index=['PMID', 'Sentence'], 
                             columns='Class', 
                             values='Entity', 
                             aggfunc=lambda x: ', '.join(x),
                             fill_value='')

# Reset index to make 'PMID' and 'Sentence' columns again
pivot_table = pivot_table.reset_index()

# Rename columns for clarity
pivot_table.columns.name = None
pivot_table.columns = ['PMID', 'Sentence', 'Activity', 'EDC', 'Target']

# Re-arrange columns to EDC, Activity, Target, PMID, Sentence
pivot_table = pivot_table[['EDC', 'Activity', 'Target', 'PMID', 'Sentence']]

# Function to split and expand the rows
def expand_rows(df):
    expanded_rows = []
    
    for _, row in df.iterrows():
        edc_list = row['EDC'].split(', ')
        target_list = row['Target'].split(', ')
        
        for edc, target in product(edc_list, target_list):
            expanded_rows.append({
                'EDC': edc,
                'Activity': row['Activity'],
                'Target': target,
                'PMID': row['PMID'],
                'Sentence': row['Sentence']
            })
    
    return pd.DataFrame(expanded_rows)

# Apply the function to expand rows
expanded_df = expand_rows(pivot_table)

# Drop duplicates based on EDC, Target, and PMID
expanded_df_no_activity = expanded_df.drop_duplicates(subset=['EDC', 'Target', 'PMID'])

# Group by EDC and Target to get counts and aggregate PMIDs
interaction_summary = expanded_df_no_activity.groupby(['EDC', 'Target']).agg(
    Count=('PMID', 'size'),
    PMIDs=('PMID', lambda x: ', '.join(map(str, x)))
).reset_index()

# Save the interaction summary to a CSV file
interaction_summary.to_csv('interaction_summary.csv', index=False)

# Drop duplicates based on EDC, Target, PMID, and Activity
expanded_df_with_activity = expanded_df.drop_duplicates(subset=['EDC', 'Target', 'PMID', 'Activity'])

# Group by EDC, Target, and Activity to get counts and aggregate PMIDs
interaction_summary_with_activity = expanded_df_with_activity.groupby(['EDC', 'Target', 'Activity']).agg(
    Count=('PMID', 'size'),
    PMIDs=('PMID', lambda x: ', '.join(map(str, x)))
).reset_index()

# Rearrange columns to place Activity between EDC and Target
interaction_summary_with_activity = interaction_summary_with_activity[['EDC', 'Activity', 'Target', 'Count', 'PMIDs']]

# Save the summary with activity to a TSV file
interaction_summary_with_activity.to_csv('interaction_summary_with_activity.tsv', sep='\t', index=False)
