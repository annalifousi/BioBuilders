import nltk
import csv
from nltk.tokenize import sent_tokenize

# Ensure NLTK data is downloaded
nltk.download('punkt')

# Parameters
percentage = 15  # Change this value to select a different percentage
csv_file_path = 'articles.csv'  # Update this path

# Create the output file name with the percentage
output_file_path = f'annotation_sentences_{percentage}percent.txt'

# Read the CSV file and extract abstracts
abstracts = []
with open(csv_file_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        abstract = row['Abstract']
        abstracts.append(abstract)

# Tokenize each abstract into sentences
sentences = []
for abstract in abstracts:
    sentences.extend(sent_tokenize(abstract))

# Calculate the specified percentage of the total number of sentences
num_sentences = len(sentences)
selected_count = max(1, num_sentences * percentage // 100)  # Ensure at least one sentence is selected

# Select the specified percentage of sentences
selected_sentences = sentences[:selected_count]

# Save selected sentences to a text file
with open(output_file_path, 'w') as outfile:
    for sentence in selected_sentences:
        outfile.write(sentence + "\n")

print(f"Selected {percentage}% of sentences and saved to {output_file_path}")
