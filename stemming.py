import pandas as pd
import nltk

# Download the necessary NLTK resources
nltk.download('punkt')


# Define a function to clean and tokenize the text
def clean_and_tokenize(text):
    if pd.isna(text):
        text = ""
    # Remove specified punctuations
    text = text.translate(str.maketrans('', '', "}{}[]!@#$%^&*()§±;<>?:_"))
    # Convert text to lowercase
    text = text.lower()
    # Tokenize the text into sentences
    sentences = nltk.sent_tokenize(text)
    return sentences


# Load the CSV file
df = pd.read_csv('articles.csv')

# Initialize an empty list to store the processed data
processed_data = []

# Process each row in the DataFrame
for index, row in df.iterrows():
    title_sentences = clean_and_tokenize(row['Title'])
    abstract_sentences = clean_and_tokenize(row['Abstract'])

    # Combine title and abstract sentences
    all_sentences = title_sentences + abstract_sentences

    # Add each sentence with its ID to the processed_data list
    for sentence_id, sentence in enumerate(all_sentences, 1):
        processed_data.append({
            'PMID': row['PMID'],
            'PMCID': row['PMC ID'],
            'Authors': row['Authors'],
            'Title': row['Title'],
            'Date': row['Date'],
            'Sentence': sentence,
            'Sentence id': sentence_id
        })

# Create a new DataFrame from the processed data
processed_df = pd.DataFrame(processed_data)

# Save the processed data to a new CSV file
processed_df.to_csv('processed_articles.csv', index=False)


print("Processed text has been saved in processed_articles.csv")
