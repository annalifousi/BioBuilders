import os
from Bio import Entrez
import openpyxl
from openpyxl import Workbook

# Identifying the connected user
Entrez.email = 'annalifousi@gmail.com'
Entrez.tool = 'Demoscript'

# Define the query
query = '("EDC"[Title/Abstract] or "endocrine disrupting chemicals"[Title/Abstract] and "receptors"[Title/Abstract])'

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

# Determine the current working directory and file path for XML
current_dir = os.getcwd()
file_path = os.path.join(current_dir, 'venv/articles_tool.xml')

# Write records in XML file
with open(file_path, 'wb') as f:
    f.write(fetch.read())

fetch.close()  # Close the fetch object to release resources

# Print the number of articles and the file path
print(f'Number of articles returned: {num_articles}')
print(f'The file has been saved to: {file_path}')

# Define the Excel file path
excel_file_path = os.path.join(current_dir, 'venv/query_results.xlsx')

# Create a new Excel workbook if it doesn't exist, otherwise load it
if not os.path.exists(excel_file_path):
    wb = Workbook()
    ws = wb.active
    ws.title = "Query Results"
    ws.append(["Query", "Number of Articles"])
else:
    wb = openpyxl.load_workbook(excel_file_path)
    ws = wb.active

# Append the new query and number of articles to the workbook
ws.append([query, num_articles])

# Save the workbook
wb.save(excel_file_path)

print(f'The Excel file has been updated and saved to: {excel_file_path}')