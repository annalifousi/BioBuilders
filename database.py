import sqlite3
import pandas as pd
from sqlalchemy import create_engine

# Define database path
db_path = 'ner_results.db'

# Connect to SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Create or update the ner_entities table
cursor.execute('''
CREATE TABLE IF NOT EXISTS ner_entities (
    Entity_EDC TEXT,
    Entity_TARGET TEXT,
    Count INTEGER,
    articles INTEGER
)
''')

# Create or update the sentences table without the catalog column
cursor.execute('''
CREATE TABLE IF NOT EXISTS sentences (
    activity TEXT,
    endocrine_disrupting_chemical TEXT,
    target TEXT,
    counts INTEGER,
    articles INTEGER,
    PRIMARY KEY (activity, endocrine_disrupting_chemical, target)
)
''')

# Commit the changes and close the connection
conn.commit()
conn.close()

# Load the CSV data
file_path = 'final_dataframe.csv'
df = pd.read_csv(file_path)

# Display the DataFrame to verify the content
print(df.head())

# Prepare data for insertion into sentences table
sentences_df = df[['ACTIVITY', 'ENDOCRINE_DISRUPTING_CHEMICAL', 'TARGET', 'counts', 'articles']]
sentences_df.columns = ['activity', 'endocrine_disrupting_chemical', 'target', 'counts', 'articles']

# Create an SQLite database connection using SQLAlchemy
engine = create_engine(f'sqlite:///{db_path}')

# Insert data into the sentences table
sentences_df.to_sql('sentences', engine, if_exists='replace', index=False)

print("Data inserted into the SQL database successfully.")
