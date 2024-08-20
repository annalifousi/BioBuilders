import sqlite3
import pandas as pd
from sqlalchemy import create_engine

# Define the path to the SQLite database
db_path = 'ner_results.db'

# Create or connect to the SQLite database
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Create tables with the correct schema
# Adjusting based on the columns in the prepared dataframes
cursor.execute('''
CREATE TABLE IF NOT EXISTS interaction_summary (
    endocrine_disrupting_chemical TEXT,
    target TEXT,
    counts INTEGER,
    pmids TEXT
)
''')

cursor.execute('''
CREATE TABLE IF NOT EXISTS interaction_summary_with_activity (
    endocrine_disrupting_chemical TEXT,
    activity TEXT,
    target TEXT,
    pmids TEXT,
    sentences TEXT
)
''')

# Commit the changes
conn.commit()
conn.close()

# Reload and parse the TSV and CSV data correctly
interaction_summary_df = pd.read_csv('interaction_summary.tsv', sep='\t')
interaction_summary_with_activity_df = pd.read_csv('interaction_summary_with_activity.tsv', sep='\t')

# Display the columns to verify correctness
print("Columns in interaction_summary DataFrame:")
print(interaction_summary_df.columns)

print("Columns in interaction_summary_with_activity DataFrame:")
print(interaction_summary_with_activity_df.columns)

# Prepare data for insertion into the tables
interaction_summary_prepared = interaction_summary_df[['EDC', 'Target', 'Count', 'PMIDs']]
interaction_summary_prepared.columns = ['endocrine_disrupting_chemical', 'target', 'counts', 'pmids']

interaction_summary_with_activity_prepared = interaction_summary_with_activity_df[['EDC', 'Activity', 'Target', 'PMIDs', 'Sentences']]
interaction_summary_with_activity_prepared.columns = [
    'endocrine_disrupting_chemical', 
    'activity', 
    'target', 
    'pmids', 
    'sentences'
]

# Create an SQLite database connection using SQLAlchemy
engine = create_engine(f'sqlite:///{db_path}')

# Insert data into the interaction_summary table
interaction_summary_prepared.to_sql('interaction_summary', engine, if_exists='replace', index=False)

# Insert data into the interaction_summary_with_activity table
interaction_summary_with_activity_prepared.to_sql('interaction_summary_with_activity', engine, if_exists='replace', index=False)

print("Data inserted into the SQL database successfully.")
