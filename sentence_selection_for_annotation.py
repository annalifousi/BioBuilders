import nltk
import csv
import pandas as pd
from nltk.tokenize import sent_tokenize, word_tokenize
import re

# Ensure NLTK data is downloaded
nltk.download('punkt')

# Increase the CSV field size limit
import sys
import csv
csv.field_size_limit(sys.maxsize)  # This sets the field size limit to the maximum allowed by the system

# Parameters
percentage = 15  # Change this value to select a different percentage
csv_file_path = 'pmc_full_texts_with_dois_and_sentences.csv'  # Update this path
combined_catalog_path = 'combined_edc_catalog.tsv'  # Update this path
chemical_targets_path = 'combined_target_catalog.tsv'  # Update this path
#additional_words_path = 'edc_vocabulary.txt'  # Path to the additional words file (optional)

# Create the output file name with the percentage
output_file_path = f'annotation_sentences_{percentage}percent.txt'

# Optional feature switch: Set to True to include additional words, False to exclude
include_additional_words = True  # Change this to False to exclude additional words

# Function to load a TSV file into a list and handle empty files
def load_tsv_file(file_path):
    try:
        data = pd.read_csv(file_path, sep='\t', header=None)
        if data.empty:
            print(f"Warning: {file_path} is empty.")
            return []
        return data.squeeze().tolist()
    except pd.errors.EmptyDataError:
        print(f"Error: {file_path} is empty or could not be read.")
        return []

# Load the catalogs
combined_catalog = load_tsv_file(combined_catalog_path)
chemical_targets = load_tsv_file(chemical_targets_path)

# Combine both catalogs and convert all words to lowercase for case-insensitive matching
catalog_words = set(word.lower() for word in combined_catalog + chemical_targets)

# Load additional words if the feature is enabled
"""if include_additional_words:
    additional_words = load_tsv_file(additional_words_path)
    catalog_words.update(word.lower() for word in additional_words)
"""
# Read the CSV file and extract abstracts and sections
data = pd.read_csv(csv_file_path)
sections_of_interest = ["Results", "Discussion"]

# Filter the DataFrame to include only rows from the "Results" and "Discussion" sections
filtered_data = data[data['Section'].isin(sections_of_interest)]

# Tokenize each abstract into sentences and filter sections
sentences = []
for _, row in filtered_data.iterrows():
    abstract = row['Sentence']
    sentences.extend(sent_tokenize(abstract))

# Function to check if a sentence contains any catalog word
def contains_catalog_word(sentence, catalog_words):
    words = word_tokenize(sentence.lower())  # Tokenize and convert to lowercase
    return any(word in catalog_words for word in words)

# Function to check if a sentence contains only numbers and symbols
def is_valid_sentence(sentence):
    # Regular expression to check if a sentence contains only numbers and symbols
    return not re.fullmatch(r'[^\w\s]+|\d+', sentence)

# Filter sentences that contain at least one catalog word and are valid
filtered_sentences = [sentence for sentence in sentences if contains_catalog_word(sentence, catalog_words) and is_valid_sentence(sentence)]

# Calculate the specified percentage of the filtered sentences
num_filtered_sentences = len(filtered_sentences)
selected_count = max(1, num_filtered_sentences * percentage // 100)  # Ensure at least one sentence is selected

# Select the specified percentage of filtered sentences
selected_sentences = filtered_sentences[:selected_count]

# Save selected sentences to a text file
with open(output_file_path, 'w') as outfile:
    for sentence in selected_sentences:
        outfile.write(sentence + "\n")

print(f"Selected {percentage}% of sentences with catalog words and saved to {output_file_path}")
