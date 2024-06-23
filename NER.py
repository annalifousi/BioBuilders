"""
import pandas as pd
import spacy
import scispacy
import nltk
nltk.download('punkt')


import en_ner_bionlp13cg_md

nlp_bi = en_ner_bionlp13cg_md.load()
import en_ner_craft_md

nlp_cr = en_ner_craft_md.load()
import en_ner_bc5cdr_md

nlp_bc = en_ner_bc5cdr_md.load()

#nlp_core = spacy.load("en_core_sci_sm")


from bs4 import BeautifulSoup


def ner(text, pmcid, table, f):
    if f == "bi":
        doc = nlp_bi(text)
    elif f == "bc":
        doc = nlp_bc(text)
    else:
        print("ERROR")
        return table

    ent = {}
    for x in doc.ents:
        ent[x.text] = x.label_
    for k in ent:
        # count +=1
        table["ID"].append(pmcid)
        table["Entity"].append(k)
        table["Class"].append(ent[k])

    return table

def segmentation(text, meta, pmcid, article_title):
    tokens = nltk.sent_tokenize(text)
    for t in tokens:
        # print(t, "\n")
        doc = nlp_core(t)
        # print(tokens)
        meta["SENTENCE"].append(t)
        meta["DOCENTS"].append(doc.ents)
        meta["PMCID"].append(pmcid)
        meta["TITLE"].append(article_title)

    return meta



with open('articles_tool.xml', "r") as xml_file:
    soup = BeautifulSoup(xml_file, 'lxml')
all = soup.findAll('article')

f = open('metadata.txt', 'w')
par_text = open('paragraphs.txt', 'w')



table = {"ID": [], "Entity": [], "Class": []}
meta = { "PMCID": [], "TITLE": [], "SENTENCE": [], "DOCENTS": []}
paragraph_list = ""


for article in all:

    #pmcid = article.find("article-id", attrs={"pub-id-type": "pmc"}).text
    article_tag = article.find("article-title")
    article_title = article_tag.text
    paragraph_list = article.findAll('p')


    segmentation(article_title, meta, pmcid, article_title)
    ner(article_title, pmcid, table, "bi")
    ner(article_title, pmcid, table, "bc")


    if paragraph_list:
        for a in paragraph_list:
            paragraph = a.text

            segmentation(paragraph, meta, pmcid, article_title)
            ner(paragraph, pmcid, table, "bi")
            ner(paragraph, pmcid, table, "bc")



trans_df = pd.DataFrame(table)
trans_df.to_csv('Entity_bi.csv', index=False)


trans_df = pd.DataFrame(meta)
trans_df.to_csv('meta_info.csv', index=False)

"""
import pandas as pd
import spacy
from scispacy.linking import EntityLinker
import en_ner_bionlp13cg_md
import en_ner_bc5cdr_md
from scispacy.abbreviation import AbbreviationDetector

# Load the models
print("Loading models...")
nlp_bi = en_ner_bionlp13cg_md.load()
nlp_bc = en_ner_bc5cdr_md.load()
print("Models loaded successfully.")

# Initialize the Entity linker
print("Initializing EntityLinker...")
# linker = EntityLinker(name="umls", k=30)

# Add EntityLinker and AbbreviationDetector to each model's pipeline
nlp_bi.add_pipe("abbreviation_detector")
nlp_bi.add_pipe("scispacy_linker", config={"resolve_abbreviations": True, "linker_name": "umls"})
nlp_bc.add_pipe("abbreviation_detector")
nlp_bc.add_pipe("scispacy_linker", config={"resolve_abbreviations": True, "linker_name": "umls"})

print("EntityLinker initialized successfully.")

# Function to perform NER, standardize abbreviations, and add results to the table
def ner(text, pmcid, table, f):
    if f == "bi":
        doc = nlp_bi(text)
    elif f == "bc":
        doc = nlp_bc(text)
    else:
        print("ERROR")
        return table

    # Process abbreviations
    abbreviations = {}
    for abrv in doc._.abbreviations:
        long_form = str(abrv._.long_form)
        short_form = str(abrv)
        text = text.replace(short_form, long_form)
        abbreviations[short_form.lower()] = long_form
    
    # Lowercase the text
    text = text.lower()
    
    # Reprocess the document with expanded abbreviations
    doc = nlp_bi(text) if f == "bi" else nlp_bc(text)

    ent = {}
    for x in doc.ents:
        ent[x.text.lower()] = x.label_
    
    for k in ent:
        table["ID"].append(pmcid)
        table["Entity"].append(k)
        table["Class"].append(ent[k])
        # Add long form if abbreviation exists, otherwise empty string
        if k in abbreviations: # long form for the chemicals when needed
            table["LongForm"].append(abbreviations[k])
        else:
            table["LongForm"].append("")

    return table

# Read the CSV file
print("Reading CSV file...")
df = pd.read_csv('articles.csv')
print("CSV file read successfully.")

# Replace NaN values with empty strings
df = df.fillna("")

df_subset = df.head(4)
# Initialize the table for storing entities
table = {"ID": [], "Entity": [], "Class": [], "LongForm": []}

# Iterate through each row in the DataFrame
for index, row in df_subset.iterrows():
    pmcid = row['PMC ID']
    title = row['Title']
    abstract = row['Abstract']

    # Ensure the text fields are strings
    title = str(title)
    abstract = str(abstract)

    # Print progress
    print(f"Processing row {index + 1}/{len(df)} - PMC ID: {pmcid}")

    # Apply NER to the title
    ner(title, pmcid, table, "bi")
    ner(title, pmcid, table, "bc")

    # Apply NER to the abstract
    ner(abstract, pmcid, table, "bi")
    ner(abstract, pmcid, table, "bc")

# Convert the results table to a DataFrame and save to a CSV file
entity_df = pd.DataFrame(table)
entity_df.to_csv('ner_data.csv', index=False)
print("Entity extraction completed and results saved to ner_data.csv.")
