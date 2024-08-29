import pandas as pd

# Paths to the files (Update these paths as necessary)
csv_file_path = 'processed_articles_full_text.csv'  # Path to your CSV file
combined_catalog_path = 'combined_edc_catalog.tsv'  # Path to combined EDC catalog
chemical_targets_path = 'combined_target_catalog.tsv'  # Path to chemical targets catalog
edc_vocabulary_path = 'edc_vocabulary.txt'  # Path to EDC vocabulary

# Load the CSV file
df_csv = pd.read_csv(csv_file_path)

# Load the TSV files containing the relevant words
df_combined = pd.read_csv(combined_catalog_path, sep='\t', header=None)
df_chemical_targets = pd.read_csv(chemical_targets_path, sep='\t', header=None)
df_edc_vocabulary = pd.read_csv(edc_vocabulary_path, sep='\t', header=None)

# Combine all relevant words into separate sets
relevant_words_combined = set(df_combined[0].tolist())
relevant_words_chemical_targets = set(df_chemical_targets[0].tolist())
relevant_words_edc_vocabulary = set(df_edc_vocabulary[0].tolist())

# Function to check if a sentence contains at least one word from each set
def contains_words_from_all(sentence, words_combined, words_chemical_targets, words_edc_vocabulary):
    sentence_lower = sentence.lower()
    contains_combined = any(word in sentence_lower for word in words_combined)
    contains_chemical_targets = any(word in sentence_lower for word in words_chemical_targets)
    contains_edc_vocabulary = any(word in sentence_lower for word in words_edc_vocabulary)
    return contains_combined and contains_chemical_targets and contains_edc_vocabulary

# Filter sentences that contain at least one word from all three sets
filtered_sentences = df_csv[df_csv['Sentence'].apply(lambda x: contains_words_from_all(x, relevant_words_combined, relevant_words_chemical_targets, relevant_words_edc_vocabulary))]

# Sample 5% of the filtered sentences
sampled_sentences = filtered_sentences.sample(frac=0.05, random_state=1)

# Extract only the 'Sentence' column and clean up any extraneous quotation marks
sentences_only = sampled_sentences['Sentence'].str.strip().str.replace('"', '', regex=False)

# Remove duplicate sentences
sentences_unique = sentences_only.drop_duplicates()

# Save the result to a tab-separated text file
output_file_path = 'filtered_sampled_sentences.txt'
with open(output_file_path, 'w') as file:
    for sentence in sentences_unique:
        file.write(f"{sentence}\n")

print(f"Number of sentences after filtering: {len(filtered_sentences)}")
print(f"Number of sentences in the sample: {len(sentences_unique)}")
print(f"Filtered, sampled, and deduplicated sentences saved to {output_file_path}")
