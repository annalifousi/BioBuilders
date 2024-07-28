import pandas as pd
import os

# Parameters
input_csv_file = 'processed_articles.csv'  # Path to your CSV file
vocabulary_file = 'edc_vocabulary.tsv'  # Path to your vocabulary file
catalog_file = 'combined_edc_catalog.tsv'  # Path to the combined EDC catalog file
receptors_file = 'receptors.tsv'  # Path to the receptors file (assumed TSV)
output_dir = 'output'  # Directory to save the .txt files
max_lines_per_file = 1000  # Maximum number of lines per file
required_sentences = 100  # Minimum number of sentences required
dummy_sentences = [
    "This is a dummy sentence for training purposes.",
    "Another example of a dummy sentence to fill the dataset.",
    "Dummy sentences help to ensure the dataset has enough content.",
    "Training models often require a variety of example sentences.",
    "Adding dummy data can be useful for balancing datasets."
]  # Example dummy sentences

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load the vocabulary words
vocabulary_df = pd.read_csv(vocabulary_file, delimiter='\t', header=None)
vocabulary_words = set(vocabulary_df[0].str.lower())  # Convert to set for faster lookup

# Load the combined EDC catalog words
catalog_df = pd.read_csv(catalog_file, delimiter='\t', header=None)
catalog_words = set(catalog_df[0].str.lower())  # Convert to set for faster lookup

# Load the receptors words
receptors_df = pd.read_csv(receptors_file, delimiter='\t', header=None)
receptors_words = set(receptors_df[0].str.lower())  # Convert to set for faster lookup

# Read the CSV file
df = pd.read_csv(input_csv_file)

# Function to write text to files
def write_to_file(text_lines, file_index):
    file_path = os.path.join(output_dir, f'annotated_passages_{file_index}.txt')
    with open(file_path, 'w') as f:
        f.write('\n'.join(text_lines))
    print(f'File saved: {file_path}')

# Create a dictionary to track sentences for each vocabulary word
word_to_sentence = {word: None for word in vocabulary_words}

# Gather sentences for each vocabulary word
valid_sentences = []
for index, row in df.iterrows():
    sentence = row['Sentence']
    words_in_sentence = set(sentence.lower().split())

    # Find intersection with vocabulary words
    vocab_intersect = words_in_sentence.intersection(vocabulary_words)
    catalog_intersect = words_in_sentence.intersection(catalog_words)
    receptors_intersect = words_in_sentence.intersection(receptors_words)

    # Ensure at least one word from each set is present in the sentence
    if vocab_intersect and catalog_intersect and receptors_intersect:
        valid_sentences.append(sentence)
        # Store the sentence for each vocabulary word covered
        for word in vocab_intersect:
            if word_to_sentence[word] is None:
                word_to_sentence[word] = sentence

# Ensure we have at least `required_sentences` total, including dummy sentences
if len(valid_sentences) < required_sentences:
    # Add dummy sentences to fill the gap
    additional_needed = required_sentences - len(valid_sentences)
    valid_sentences.extend(dummy_sentences[:additional_needed])

# Shuffle the list of sentences to mix valid and dummy sentences
import random
random.shuffle(valid_sentences)

# Write sentences to files, ensuring the maximum number per file
text_lines = []
file_index = 1

for sentence in valid_sentences:
    text_lines.append(sentence)

    # Check if the number of lines has reached the maximum for a file
    if len(text_lines) >= max_lines_per_file:
        write_to_file(text_lines, file_index)
        text_lines = []
        file_index += 1

# Write any remaining lines to a final file
if text_lines:
    write_to_file(text_lines, file_index)

# Find and print words that couldn't be found and those that were found
missing_words = [word for word, sentence in word_to_sentence.items() if sentence is None]
found_words = [word for word, sentence in word_to_sentence.items() if sentence is not None]

print(f"Total vocabulary words: {len(vocabulary_words)}")
print(f"Words found: {len(found_words)}")
print(f"Words not found: {len(missing_words)}")

if missing_words:
    print(f"Words not found in any sentence: {', '.join(missing_words)}")
else:
    print("All words were found in the sentences.")

print("Processing complete.")
