import csv
import re


def extract_entities(sentence, tagged_words):
    entities = []
    for tagged_word in tagged_words:
        if '(' in tagged_word and ')' in tagged_word:
            word, tag = tagged_word.split('(')
            tag = tag.rstrip(')')
            start = sentence.find(word)
            if start != -1:
                end = start + len(word)
                entities.append((start, end, tag))
    return entities


def create_spacy_training_data(csv_file_path):
    TRAIN_DATA = []

    with open(csv_file_path, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            sentence = row['Sentence']
            tagged_words = row['Common Words in Vocabulary 1'].split(', ') + row['Common Words in Vocabulary 2'].split(
                ', ')
            entities = extract_entities(sentence, tagged_words)
            TRAIN_DATA.append((sentence, {'entities': entities}))

    return TRAIN_DATA


def save_training_data_to_py_file(training_data, output_file_path):
    with open(output_file_path, 'w', encoding='utf-8') as file:
        file.write('TRAIN_DATA = [\n')
        for sentence, annotation in training_data:
            file.write(f"    ({repr(sentence)}, {annotation}),\n")
        file.write(']\n')


if __name__ == '__main__':
    csv_file_path = 'report.csv'  # Replace with your file path
    output_file_path = 'spacy_training_data.py'  # The output file path
    training_data = create_spacy_training_data(csv_file_path)

    # Save the training data to a Python file
    save_training_data_to_py_file(training_data, output_file_path)

    # Print a message to indicate completion
    print(f'Training data saved to {output_file_path}')

