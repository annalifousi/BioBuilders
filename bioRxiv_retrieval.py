import pandas as pd
import requests


def fetch_paper_by_doi(doi, response_format="json"):
    # Construct the correct API URL
    url = f"https://api.biorxiv.org/details/biorxiv/[DOI]/na/json"

    # Send GET request to bioRxiv API
    response = requests.get(url)

    # Check the status code
    if response.status_code == 200:
        # Print the content type to verify the format of the response
        content_type = response.headers.get("Content-Type", "")
        print(f"Content-Type: {content_type}")

        if "json" in content_type:
            try:
                json_data = response.json()
                print("JSON Response:", json_data)  # Print the entire JSON response

                if 'collection' in json_data and json_data['collection']:
                    paper_info = json_data['collection'][0]
                    title = paper_info.get('title', 'No title found')
                    abstract = paper_info.get('abstract', 'No abstract found')
                    print(f"Title: {title}\nAbstract: {abstract}\n")
                else:
                    print(f"No data found for DOI {doi}.")
            except (KeyError, IndexError) as e:
                print(f"Error extracting data from JSON response: {e}")
        else:
            print(f"Unexpected content type: {content_type}")
            print(f"Response content: {response.text}")
    else:
        print(f"Error fetching data for DOI {doi}: {response.status_code}, {response.text}")


# Read the CSV file containing DOIs
df = pd.read_csv("pmc_full_texts_with_dois_and_sentences.csv")

# Strip any extra spaces from column names
df.columns = df.columns.str.strip()

# Iterate over the DOIs and fetch paper details
for doi in df['DOI']:
    fetch_paper_by_doi(doi)
