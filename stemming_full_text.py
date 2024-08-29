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
receptors = pd.read_csv('combined_target_catalog.tsv', sep='\t')
targets = pd.read_csv('standardized_edc_targets.tsv', sep='\t')

# Combine the entity names from all files into a set for fast membership checking
entities = set(combined_edc_catalog['Name'].str.lower()).union(receptors['Name'].str.lower()).union(
    set(targets['Name'].str.lower()))

# Load the articles CSV file with the new format
df = pd.read_csv('pmc_sections.csv')

# Initialize lists to store the processed data and sentences without punctuation
processed_data = []
sentences_no_punctuation = []

# Initialize variables to track the current article and sentence ID
current_pmc_id = None
sentence_id_counter = 1
combined_text = ""

# Process each row in the DataFrame
for index, row in df.iterrows():
    # Check if we are processing a new article based on PMC ID
    if row['PMC ID'] != current_pmc_id:
        # If it's a new article, reset the sentence ID counter and store combined text
        if current_pmc_id is not None:
            # Store combined text for the previous article
            sentences_no_punctuation.append({
                'PMC ID': current_pmc_id,
                'DOI': current_doi,
                'Authors': current_authors,
                'Title': current_title,
                'Publication Date': current_pub_date,
                'Text': combined_text
            })

        # Update current article metadata
        current_pmc_id = row['PMC ID']
        current_doi = row['DOI']
        current_authors = row['Authors']
        current_title = row['Title']
        current_pub_date = row['Publication Date']

        # Reset sentence counter and combined text
        sentence_id_counter = 1
        combined_text = ""

    # Clean and tokenize the article section text
    section_sentences, section_no_punctuation = clean_and_tokenize(row['Article'])

    # Append the section text to the combined text for the article
    combined_text += " " + section_no_punctuation

    # Filter sentences containing any entity from combined_edc_catalog or receptors
    filtered_sentences = [sentence for sentence in section_sentences if any(entity in sentence for entity in entities)]

    # Add each filtered sentence with its ID to the processed_data list
    for sentence in filtered_sentences:
        processed_data.append({
            'PMC ID': row['PMC ID'],
            'DOI': row['DOI'],
            'Authors': row['Authors'],
            'Title': row['Title'],
            'Publication Date': row['Publication Date'],
            'Section': row['Section'],
            'Sentence': sentence,
            'Sentence id': sentence_id_counter
        })
        sentence_id_counter += 1

# Add the last article's combined text to sentences_no_punctuation
if current_pmc_id is not None:
    sentences_no_punctuation.append({
        'PMC ID': current_pmc_id,
        'DOI': current_doi,
        'Authors': current_authors,
        'Title': current_title,
        'Publication Date': current_pub_date,
        'Text': combined_text
    })

# Create DataFrames from the processed data and sentences without punctuation
processed_df = pd.DataFrame(processed_data)
sentences_no_punctuation_df = pd.DataFrame(sentences_no_punctuation)

# Save the processed data to new CSV files
processed_df.to_csv('processed_articles_full_text.csv', index=False)
sentences_no_punctuation_df.to_csv('sentences_no_punctuation.csv', index=False)

print("Processed text has been saved in processed_articles.csv")
print("Sentences without punctuation have been saved in sentences_no_punctuation.csv")
