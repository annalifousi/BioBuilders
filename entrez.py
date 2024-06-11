import os
from Bio import Entrez
import openpyxl
from openpyxl import Workbook

# Identifying the connected user
Entrez.email = 'annalifousi@gmail.com'
Entrez.tool = 'Demoscript'

# Define the query
#query = '("EDC"[Title/Abstract] or "endocrine disrupting chemicals"[Title/Abstract] and "receptors"[Title/Abstract])'
#query = '("endocrine disrupting chemicals"[Title/Abstract] and "receptors"[Title/Abstract])'
#query = '("EDC"[Title/Abstract] or "endocrine disrupting chemicals"[Title/Abstract] and "receptors"[Title/Abstract] or "receptor"[Title/Abstract])'
#query = '("EDC"[Title/Abstract] or "endocrine disrupting chemicals"[Title/Abstract] and "receptors"[Title/Abstract] or "receptor"[Title/Abstract])'
#query = '("Endocrine Disruptors"[MeSH Terms] OR "endocrine disrupting chemicals"[Title/Abstract] AND "receptors"[Title/Abstract] OR "receptor"[Title/Abstract]) NOT review[Publication Type]'
query = ('("Endocrine Disruptors"[MeSH Terms] OR "endocrine disrupting chemicals"[Title/Abstract] OR "EDC"[Title/Abstract]) '
         'AND ("Receptors, Endocrine"[MeSH Terms] OR "receptors"[Title/Abstract] OR "receptor"[Title/Abstract]) '
         'AND ("binding"[Title/Abstract] OR "interaction"[Title/Abstract] OR "affinity"[Title/Abstract] OR "assay"[Title/Abstract] OR "experiment"[Title/Abstract]) '
         'NOT review[Publication Type]')


# Searching for the query in Entrez
info = Entrez.esearch(db="pubmed", retmax=100000, term=query)
# Parsing the XML data
record = Entrez.read(info)

# Count the number of articles
num_articles = len(record['IdList'])

# Retrieve records in XML format
fetch = Entrez.efetch(db='pubmed',
                      retmode='xml',
                      id=record['IdList'],
                      rettype='full')

# Determine the file path for XML
file_path = '/Users/annalifousihotmailcom/python/DTU Biobuilders/venv/articles_tool.xml'

# Write records in XML file
with open(file_path, 'wb') as f:
    f.write(fetch.read())

fetch.close()  # Close the fetch object to release resources

# Print the number of articles and the file path
print(f'Number of articles returned: {num_articles}')
print(f'The file has been saved to: {file_path}')

# Define the Excel file path
excel_file_path = '/Users/annalifousihotmailcom/python/DTU Biobuilders/venv/query_results.xlsx'

# Create a new Excel workbook if it doesn't exist, otherwise load it
if not os.path.exists(excel_file_path):
    wb = Workbook()
    ws = wb.active
    ws.title = "Query Results"
    ws.append(["Query", "Number of Articles"])
else:
    wb = openpyxl.load_workbook(excel_file_path)
    ws = wb.active

# Check the last entry to avoid duplication
last_row = ws.max_row
if last_row > 1:  # Ensure there's more than just the header
    last_query = ws.cell(row=last_row, column=1).value
    last_num_articles = ws.cell(row=last_row, column=2).value
    if last_query == query and last_num_articles == num_articles:
        print("The query and number of articles are the same as the last entry. No new entry added.")
    else:
        ws.append([query, num_articles])
        print(f'New entry added: {query} - {num_articles}')
else:
    ws.append([query, num_articles])
    print(f'New entry added: {query} - {num_articles}')

# Save the workbook
wb.save(excel_file_path)

print(f'The Excel file has been updated and saved to: {excel_file_path}')
