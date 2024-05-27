### Takes too long -- hundreds of thousands of reviewed results, millions total for each chebi id

import requests
import re
from requests.adapters import HTTPAdapter, Retry
import pandas as pd


# Initializing retry configuration for HTTP requests
retries = Retry(
    total=5,
    backoff_factor=0.25,
    status_forcelist=[500, 502, 503, 504]
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# Initializing pagination (used for gathering ec numbers from Uniprot)
def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

        
re_next_link = re.compile(r'<([^>]+)>; rel="next"')

# Initializing retry configuration for HTTP requests
retries = Retry(
    total=5,
    backoff_factor=0.25,
    status_forcelist=[500, 502, 503, 504]
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.search(headers["Link"])
        if match:
            return match.group(1)
    return None

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = int(response.headers.get("x-total-results", 0))
        yield response.text, total
        batch_url = get_next_link(response.headers)

# Main function to fetch and process data with pagination
chebi_list = ['17234', '31440']

def chebi2ec():
    chebi2ec_df = []

    # REST API base URL
    base_url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Cec%2Corganism_name%2Corganism_id%2Ccc_cofactor&format=tsv&size=500'

    for chebi_id in chebi_list:
        url = f'{base_url}&query=%28%28taxonomy_id%3A2%29+OR+%28taxonomy_id%3A2157%29+AND+%28chebi%3A{chebi_id}%29%29+AND+%28reviewed%3Atrue%29'
        
        for batch_data, total in get_batch(url):
            lines = batch_data.splitlines()
            
            # Iterate through lines to extract EC numbers
            for line in lines[1:]:  # Skip header line
                columns = line.split('\t')
                if len(columns) > 1:
                    ec_number = columns[1]  # Assuming EC number is the second column
                    chebi2ec_df.append({"ChEBI ID": chebi_id, "Enzyme": ec_number})

    # Convert the list of dictionaries into a DataFrame
    chebi2ec_df = pd.DataFrame(chebi2ec_df)
    
    return chebi2ec_df

chebi2ec_df = chebi2ec()
print(chebi2ec_df)

### Returns all entries associated with the ChEBI identifier molecule -- can be substrate, product, etc...