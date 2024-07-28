import pandas as pd
import json
import os
from tqdm import tqdm
import spacy
from spacy.tokens import DocBin

# Load the data
with open('annotations.json', 'r') as f:
    data = json.load(f)

entity_name = "ACTIVITY"

train_data = data['annotations']
train_data = [tuple(i) for i in train_data]

for i in train_data:
    if i[1]['entities'] == []:
        i[1]['entities'] = [(0, 0, entity_name)]  # List with one tuple
    else:
        i[1]['entities'] = [tuple(entity) for entity in i[1]['entities']]

# Process the data with spaCy
nlp = spacy.load("en_core_web_sm")
db = DocBin()

skipped_count = 0

for text, annot in tqdm(train_data):
    doc = nlp.make_doc(text)
    ents = []
    for entity in annot.get("entities", []):
        if isinstance(entity, (list, tuple)) and len(entity) == 3:
            start, end, label = entity
            if start < end:  # Ensure valid span
                span = doc.char_span(start, end, label=label, alignment_mode="contract")
                if span is not None:
                    ents.append(span)
                else:
                    print(f"Skipping entity with invalid span: {entity}")
                    skipped_count += 1
            else:
                print(f"Skipping entity with zero or invalid range: {entity}")
                skipped_count += 1
        else:
            print(f"Unexpected entity format: {entity}")
            skipped_count += 1
    doc.ents = ents
    db.add(doc)

# Save the DocBin object to disk
save_path = '/Users/annalifousihotmailcom/python/BioBuilders/train.spacy'
db.to_disk(save_path)
print(f"DocBin saved to {save_path}")
print(f"Total skipped entities: {skipped_count}")
