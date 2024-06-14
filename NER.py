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