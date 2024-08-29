from Bio import Entrez
from bs4 import BeautifulSoup
import csv
import numpy as np

# Set up your email for NCBI Entrez (required)
Entrez.email = "annalifousi@gmail.com"

# List of major sections we want to extract, excluding "methods"
major_sections = ['introduction', 'results', 'discussion', 'conclusion']

# Function to fetch full text from PMC using PMC IDs
def fetch_pmc_data(pmc_ids):
    # Initialize counters
    successful_retrievals = 0
    retrieval_errors = 0
    no_full_text_count = 0

    # Prepare CSV file to save the results
    with open('pmc_sections.csv', 'w', newline='', encoding='utf-8') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['PMC ID', 'DOI', 'Title', 'Authors', 'Publication Date', 'Section', 'Article'])

        for pmc_id in pmc_ids:
            try:
                # Fetch the article from PMC
                handle = Entrez.efetch(db="pmc", id=pmc_id, retmode="xml")
                xml_content = handle.read()
                handle.close()

                # Parse XML content
                soup = BeautifulSoup(xml_content, 'xml')

                # Extract article details
                title = soup.find('article-title').text if soup.find('article-title') else 'N/A'

                authors = []
                for author in soup.find_all('contrib', {'contrib-type': 'author'}):
                    last_name = author.find('surname')
                    first_name = author.find('given-names')
                    if last_name and first_name:
                        authors.append(f"{first_name.text} {last_name.text}")
                authors_text = ', '.join(authors) if authors else 'N/A'

                pub_date = soup.find('pub-date')
                pub_date_text = pub_date.find('year').text if pub_date and pub_date.find('year') else 'N/A'

                # Extract DOI
                doi = soup.find('article-id', {'pub-id-type': 'doi'})
                doi_text = doi.text if doi else 'N/A'

                # Extract the body text and major sections
                body = soup.find('body')

                if body:
                    full_text_found = False

                    # Traverse each section within the body
                    for sec in body.find_all('sec'):
                        sec_title = sec.find('title')
                        section_name = sec_title.text.lower() if sec_title else "N/A"

                        # Collect and store the text only if it's within a recognized major section
                        if any(keyword in section_name for keyword in major_sections) and "methods" not in section_name:
                            section_text = ' '.join([p.text for p in sec.find_all('p')])
                            full_text_found = True

                            # Write the section as a separate record in the CSV
                            record = [pmc_id, doi_text, title, authors_text, pub_date_text, section_name, section_text]
                            csvwriter.writerow(record)
                            print(record)

                    if not full_text_found:
                        no_full_text_count += 1

                else:
                    # If no body found, write a row with NaN in the "Section" and "Article" columns
                    no_full_text_count += 1
                    record = [pmc_id, doi_text, title, authors_text, pub_date_text, np.nan, np.nan]
                    csvwriter.writerow(record)
                    print(f"No body text for PMC ID: {pmc_id}. Recorded as NaN.")

                print(f"Processed PMC ID: {pmc_id}")

                # Increment successful retrieval counter
                successful_retrievals += 1

            except Exception as e:
                print(f"Error retrieving PMC ID {pmc_id}: {e}")
                # Increment error counter
                retrieval_errors += 1

    # Print final counts
    print(f"Total successful retrievals: {successful_retrievals}")
    print(f"Total retrieval errors: {retrieval_errors}")
    print(f"Total articles without full text: {no_full_text_count}")

# Step 1: Read PMC IDs from the CSV
pmc_ids = []
with open('articles.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if there is one
    for row in reader:
        pmc_id = row[-1]  # Assuming the last column has the PMC ID
        if pmc_id:
            pmc_ids.append(pmc_id)

# Step 2: Fetch and save data
fetch_pmc_data(pmc_ids)
