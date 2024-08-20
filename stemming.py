import pandas as pd
import nltk

# Download the necessary NLTK resources
nltk.download('punkt')

# Define a function to clean punctuation from text
def remove_punctuation(text):
    if pd.isna(text):
        text = ""
    else:
        text = str(text)  # Ensure the text is a string
    
    # Define the punctuations you want to remove, including double quotes
    remove_punctuations = '"\'(){}[]!@#$%^&*()§±;<>?:_=®/'
    
    # Create a translation table that maps the specified punctuation characters to None
    translator = str.maketrans('', '', remove_punctuations)
    
    # Remove specified punctuations
    text_no_punctuation = text.translate(translator)
    
    # Convert text to lowercase
    text_no_punctuation = text_no_punctuation.lower()
    
    return text_no_punctuation

# Define a function to clean and tokenize the text
def clean_and_tokenize(text):
    if pd.isna(text):
        text = ""
    else:
        text = str(text)  # Ensure the text is a string
    
    # Remove punctuation
    text_no_punctuation = remove_punctuation(text)
    
    # Tokenize the text into sentences
    sentences = nltk.sent_tokenize(text_no_punctuation)
    
    return sentences, text_no_punctuation

# Load the entity names from combined_edc_catalog.tsv and receptors.tsv
combined_edc_catalog = pd.read_csv('combined_edc_catalog.tsv', sep='\t')
receptors = pd.read_csv('receptors.tsv', sep='\t')
targets = pd.read_csv('standardized_edc_targets.tsv', sep='\t')
# Combine the entity names from both files into a set for fast membership checking
entities = set(combined_edc_catalog['Name'].str.lower()).union(receptors['Name'].str.lower()).union(set(targets['Name'].str.lower()))

# Load the articles CSV file
df = pd.read_csv('articles.csv')

# Ensure that all relevant text fields are strings
df['Title'] = df['Title'].apply(lambda x: remove_punctuation(str(x)) if pd.notna(x) else "")
df['Abstract'] = df['Abstract'].apply(lambda x: remove_punctuation(str(x)) if pd.notna(x) else "")
df['Authors'] = df['Authors'].apply(lambda x: remove_punctuation(str(x)) if pd.notna(x) else "")
df['Date'] = df['Date'].apply(lambda x: remove_punctuation(str(x)) if pd.notna(x) else "")
df['PMID'] = df['PMID'].apply(lambda x: remove_punctuation(str(x)) if pd.notna(x) else "")
df['PMC ID'] = df['PMC ID'].apply(lambda x: remove_punctuation(str(x)) if pd.notna(x) else "")

# Initialize lists to store the processed data and sentences without punctuation
processed_data = []
sentences_no_punctuation = []

# Process each row in the DataFrame
for index, row in df.iterrows():
    # Clean and tokenize both title and abstract
    title_sentences, title_no_punctuation = clean_and_tokenize(row['Title'])
    abstract_sentences, abstract_no_punctuation = clean_and_tokenize(row['Abstract'])

    # Store sentences without punctuation
    sentences_no_punctuation.append({
        'PMID': row['PMID'],
        'PMCID': row['PMC ID'],
        'Authors': row['Authors'],
        'Title': row['Title'],
        'Date': row['Date'],
        'Text': title_no_punctuation + " " + abstract_no_punctuation
    })

    # Combine title and abstract sentences
    all_sentences = title_sentences + abstract_sentences

    # Filter sentences containing any entity from combined_edc_catalog or receptors
    filtered_sentences = [sentence for sentence in all_sentences if any(entity in sentence for entity in entities)]

    # Add each filtered sentence with its ID to the processed_data list
    for sentence_id, sentence in enumerate(filtered_sentences, 1):
        processed_data.append({
            'PMID': row['PMID'],
            'PMCID': row['PMC ID'],
            'Authors': row['Authors'],
            'Title': row['Title'],
            'Date': row['Date'],
            'Sentence': sentence,
            'Sentence id': sentence_id
        })

# Create DataFrames from the processed data and sentences without punctuation
processed_df = pd.DataFrame(processed_data)
sentences_no_punctuation_df = pd.DataFrame(sentences_no_punctuation)

# Save the processed data to new CSV files
processed_df.to_csv('processed_articles.csv', index=False)
sentences_no_punctuation_df.to_csv('sentences_no_punctuation.csv', index=False)

print("Processed text has been saved in processed_articles.csv")
print("Sentences without punctuation have been saved in sentences_no_punctuation.csv")