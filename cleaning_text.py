import csv
import re

def clean_text(text):
    # Convert to lowercase and remove symbols except for alphanumeric characters and spaces
    text = text.lower()
    text = re.sub(r'[^\w\s]', '', text)
    return text.strip()

def process_sentence(sentence, vocabulary1, vocabulary2):
    words = sentence.split()
    common_words1 = [(word, 'ACTIVITY') for word in words if word in vocabulary1]
    common_words2 = [(word, 'RECEPTOR') for word in words if word in vocabulary2]
    return common_words1, common_words2

def main(csv_file, vocabulary_file1, vocabulary_file2, output_file):
    # Load vocabularies
    with open(vocabulary_file1, 'r', encoding='utf-8') as f1, open(vocabulary_file2, 'r', encoding='utf-8') as f2:
        vocabulary1 = set(word.strip() for word in f1.readlines())
        vocabulary2 = set(word.strip() for word in f2.readlines())

    results = []
    articles_with_common_words = 0

    # Process CSV file
    with open(csv_file, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        processed_articles = set()  # To keep track of processed articles
        for row in reader:
            pmid = row['PMID']
            title = row['Title']
            abstract = row['Abstract']

            if pmid in processed_articles:
                continue

            # Clean and tokenize title and abstract
            title_cleaned = clean_text(title)
            abstract_cleaned = clean_text(abstract)

            # Process sentences in abstract
            sentences = re.split(r'[.!?]', abstract_cleaned)
            found_common_words = False
            for sentence in sentences:
                if sentence.strip():  # Ensure it's not an empty sentence
                    common_words1, common_words2 = process_sentence(sentence, vocabulary1, vocabulary2)
                    if common_words1 and common_words2:
                        found_common_words = True
                        results.append({
                            'PMID': pmid,
                            'Title': title_cleaned,
                            'Abstract': abstract_cleaned,
                            'Sentence': sentence,
                            'Common Words in Vocabulary 1': ', '.join([f"{word}({tag})" for word, tag in common_words1]),
                            'Common Words in Vocabulary 2': ', '.join([f"{word}({tag})" for word, tag in common_words2])
                        })

            if found_common_words:
                articles_with_common_words += 1
                processed_articles.add(pmid)

    # Write results to CSV report file
    fieldnames = ['PMID', 'Title', 'Abstract', 'Sentence', 'Common Words in Vocabulary 1', 'Common Words in Vocabulary 2']
    with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    # Write the report summary to a text file
    with open('report_summary.txt', 'w', encoding='utf-8') as summary_file:
        summary_file.write(f"Total articles with common words in both vocabularies: {articles_with_common_words}\n")

    # Print the count of articles with common words in both vocabularies
    print(f"Total articles with common words in both vocabularies: {articles_with_common_words}")

    # Return count of articles with common words in both vocabularies
    return articles_with_common_words

if __name__ == "__main__":
    csv_file = 'articles.csv'
    vocabulary_file1 = 'edc_vocabulary.txt'
    vocabulary_file2 = 'receptors.txt'
    output_file = 'report.csv'
    articles_count = main(csv_file, vocabulary_file1, vocabulary_file2, output_file)
