TRAIN_DATA = [
    ("The chemical can enhance the hormone activity.", {"entities": [(19, 27, "TARGET")]}),
    ("It is known to inhibit the receptor function.", {"entities": [(13, 20, "TARGET")]}),
    # Add more sentences with labeled entities
]
import spacy
from spacy.training.example import Example

# Load a pre-existing spacy model
nlp = spacy.blank("en")

# Create an NER pipeline component and add it to the pipeline
if "ner" not in nlp.pipe_names:
    ner = nlp.add_pipe("ner")
else:
    ner = nlp.get_pipe("ner")

# Add new labels to the NER component
for _, annotations in TRAIN_DATA:
    for ent in annotations.get("entities"):
        ner.add_label(ent[2])

# Disable other pipeline components during training
pipe_exceptions = ["ner"]
unaffected_pipes = [pipe for pipe in nlp.pipe_names if pipe not in pipe_exceptions]

# Training the model
import random
from spacy.util import minibatch, compounding

with nlp.disable_pipes(*unaffected_pipes):
    optimizer = nlp.begin_training()
    for i in range(20):
        random.shuffle(TRAIN_DATA)
        losses = {}
        batches = minibatch(TRAIN_DATA, size=compounding(4.0, 32.0, 1.001))
        for batch in batches:
            for text, annotations in batch:
                doc = nlp.make_doc(text)
                example = Example.from_dict(doc, annotations)
                nlp.update([example], drop=0.5, losses=losses)
        print("Losses", losses)

# Save the model
nlp.to_disk("edc_ner_model")
nlp = spacy.load("edc_ner_model")

# Test the model
test_text = "The compound is known to stimulate the protein synthesis."
doc = nlp(test_text)

for ent in doc.ents:
    print(ent.text, ent.label_)
