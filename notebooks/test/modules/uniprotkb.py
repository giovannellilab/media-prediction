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

# Compiling regex for extracting next link from headers
re_next_link = re.compile(r'<([^>]+)>; rel="next"')

# Function to get the next link for pagination
def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.search(headers["Link"])
        if match:
            return match.group(1)
    return None

# Initialize an empty list to store DataFrame rows
def taxa2ec(id_list: list):
    
    taxa2ec_df = []
    no_data = 0
    
    # REST API base URL
    base_url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Cec%2Corganism_name%2Corganism_id%2Ccc_cofactor%2Cid&format=tsv&size=500'

    for taxonomy_id in id_list:
        #url = f'{base_url}&query=%28organism_name%3A{taxonomy_id}%29+AND+%28ec%3A*%29' #Search by NCBI ID
        url = f'{base_url}&query=%28organism_name%3A{taxonomy_id}%29+AND+%28ec%3A*%29+AND+%28reviewed%3Atrue%29' #Search by reviewed species

        try:
            response = session.get(url)
            response.raise_for_status()
            lines = response.text.splitlines()

            # Check if no results found
            if len(lines) <= 1:
                no_data += 1
#                print(f"No data found for {taxonomy_id}")
                continue

            # Iterate through lines to extract EC numbers
            for line in lines[1:]:  # Skip header line
                columns = line.split('\t')
                if len(columns) > 1:
                    ec_number = columns[1]  # Assuming EC number is the second column
                    taxa2ec_df.append({"Taxa ID": taxonomy_id, "Enzyme": ec_number})

        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred for {taxonomy_id}: {http_err}")
        except Exception as err:
            print(f"An error occurred for {taxonomy_id}: {err}")

    # Convert the list of dictionaries into a DataFrame
    taxa2ec_df = pd.DataFrame(taxa2ec_df)
    
    print(f"{no_data} species with no data")
    
    return taxa2ec_df

def ec_info(id_list: list):
    
    info_df = []
    
    # REST API base URL
    base_url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Cec%2Corganism_name%2Corganism_id%2Ccc_cofactor%2Cid&format=tsv&size=500'

    for taxonomy_id in id_list:
        #url = f'{base_url}&query=%28organism_name%3A{taxonomy_id}%29+AND+%28ec%3A*%29' #Search by NCBI ID
        url = f'{base_url}&query=%28organism_name%3A{taxonomy_id}%29+AND+%28ec%3A*%29+AND+%28reviewed%3Atrue%29' #Search by reviewed species

        try:
            response = session.get(url)
            response.raise_for_status()
            lines = response.text.splitlines()

            # Check if no results found
            if len(lines) <= 1:
                no_data += 1
#                print(f"No data found for {taxonomy_id}")
                continue

            # Iterate through lines to extract EC numbers
            for line in lines[1:]:  # Skip header line
                columns = line.split('\t')
                if len(columns) > 1:
                    ec_number = columns[1]  # Assuming EC number is the second column
                    info_df.append({"Taxa ID": taxonomy_id, "Enzyme": ec_number})

        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred for {taxonomy_id}: {http_err}")
        except Exception as err:
            print(f"An error occurred for {taxonomy_id}: {err}")

    # Convert the list of dictionaries into a DataFrame
    info_df = pd.DataFrame(info_df)
    
    return taxa2ec_df