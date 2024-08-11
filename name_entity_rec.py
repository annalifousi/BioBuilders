import pandas as pd
import spacy
from collections import defaultdict
from scispacy.abbreviation import AbbreviationDetector
from collections import Counter
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md

# Load models
print("Loading models...")
nlp_bc = en_ner_bc5cdr_md.load()
nlp_custom = spacy.load('./model-best')
nlp_edc_custom = spacy.load('./EDC_target_model/output/EDC_model-best')

# Add AbbreviationDetector to each model's pipeline
for nlp_model in [nlp_bc, nlp_custom, nlp_edc_custom]:
    if "abbreviation_detector" not in nlp_model.pipe_names:
        nlp_model.add_pipe("abbreviation_detector")

# Load EDC and receptor terms
def load_terms(filename):
    df = pd.read_csv(filename, sep='\t')
    return set(df['Name'].str.lower())

edc_terms = load_terms('combined_edc_catalog.tsv')
receptor_terms = load_terms('standardized_edc_targets.tsv')
additional_receptor_terms = load_terms('receptors.tsv')
receptor_terms.update(additional_receptor_terms)

# Define N-gram and corpus frequency functions
def extract_ngrams(text, n=4):
    ngrams = set()
    for token in text.split():
        token = token.lower()
        for i in range(len(token) - n + 1):
            ngram = token[i:i+n]
            ngrams.add(ngram)
    return ngrams

def compute_corpus_frequency(texts):
    word_counter = Counter()
    for text in texts:
        words = text.split()
        word_counter.update(words)
    return word_counter

# Preprocess the text data
def preprocess_texts(df):
    sentences = df.groupby(['PMID', 'Sentence id'])['Sentence'].first().to_dict()
    unique_titles = df['Title'].drop_duplicates().tolist()
    unique_sentences = df['Sentence'].drop_duplicates().tolist()
    texts = unique_titles + unique_sentences
    corpus_counter = compute_corpus_frequency(texts)
    return corpus_counter, sentences

# Load the dataset
df = pd.read_csv('processed_articles.csv')

# Preprocess the texts to extract corpus frequencies
corpus_counter, sentences = preprocess_texts(df)

# Extract and save n-grams and corpus frequency terms
def extract_features(df, corpus_counter, theta_freq=10):
    feature_map = defaultdict(dict)
    for index, row in df.iterrows():
        pmid = row['PMID']
        sentence_id = row['Sentence id']
        sentence = row['Sentence']

        print(f"Processing sentence: {sentence}")  # Print the sentence being processed

        # Extract corpus frequency features
        words = sentence.split()
        for word in words:
            if corpus_counter[word.lower()] >= theta_freq:
                feature_key = f"corpus_freq_{word.lower()}"
                feature_map[(pmid, sentence_id, feature_key)] = {"Label": "UNKNOWN"}

    return pd.DataFrame([
        (pmid, sentence_id, entity_key, details["Label"])
        for (pmid, sentence_id, entity_key), details in feature_map.items()
        if "corpus_freq_" in entity_key  # Filter only for corpus frequency terms
    ], columns=["PMID", "SentenceID", "Entity", "Label"])

# Extract features from the dataset
features_df = extract_features(df, corpus_counter)

# Remove 'corpus_freq_' from the entities
features_df['Entity'] = features_df['Entity'].str.replace('corpus_freq_', '', regex=False)

# Perform POS tagging and label verbs
def label_verbs(features_df, sentences):
    annotated_entities = []
    for _, row in features_df.iterrows():
        pmid = row['PMID']
        sentence_id = row['SentenceID']
        entity = row['Entity']

        # Find the sentence that contains the entity
        sentence = sentences.get((pmid, sentence_id), "")
        doc = nlp_custom(sentence)  # Process the full sentence

        # Check if entity is in the sentence
        if any(token.text == entity for token in doc):
            pos_tagged = False
            for token in doc:
                if token.text == entity and token.pos_ == 'VERB':
                    annotated_entities.append({
                        'PMID': pmid,
                        'SentenceID': sentence_id,
                        'Entity': entity,
                        'Label': 'VERB'
                    })
                    pos_tagged = True
                    break

            if not pos_tagged:
                annotated_entities.append({
                    'PMID': pmid,
                    'SentenceID': sentence_id,
                    'Entity': entity,
                    'Label': row['Label']
                })
        else:
            annotated_entities.append({
                'PMID': pmid,
                'SentenceID': sentence_id,
                'Entity': entity,
                'Label': row['Label']
            })

    return pd.DataFrame(annotated_entities)

# Label verbs in the entities
features_df = label_verbs(features_df, sentences)

# Apply models and prioritize labels
def apply_models(features_df, models):
    annotated_entities = []
    for _, row in features_df.iterrows():
        pmid = row['PMID']
        sentence_id = row['SentenceID']
        entity = row['Entity']

        print(f"Applying models to entity: {entity}")  # Print the entity being processed by the models

        entity_lower = entity.lower()
        label = row['Label']

        # Only apply models if the label is still 'UNKNOWN'
        if label == 'UNKNOWN':
            model_labels = defaultdict(lambda: 'UNKNOWN')
            for model in models:
                doc = model(entity)
                for ent in doc.ents:
                    if ent.text.lower() == entity_lower:
                        model_labels[ent.label_] = ent.label_

            # Determine the most specific label
            final_label = 'UNKNOWN'
            if model_labels:
                final_label = max(model_labels, key=lambda k: model_labels[k])

            annotated_entities.append({
                'PMID': pmid,
                'SentenceID': sentence_id,
                'Entity': entity,
                'Label': final_label
            })
        else:
            annotated_entities.append({
                'PMID': pmid,
                'SentenceID': sentence_id,
                'Entity': entity,
                'Label': label
            })

    return pd.DataFrame(annotated_entities)

# Apply models and prioritize labels
final_df = apply_models(features_df, [nlp_bc, nlp_custom, nlp_edc_custom])

# Save the DataFrames
features_df.to_csv('features_df.csv', index=False)
features_df.to_csv('features_df_modified.csv', index=False)
final_df.to_csv('final_annotated_entities.csv', index=False)
print("Final annotated entities saved to final_annotated_entities.csv")
