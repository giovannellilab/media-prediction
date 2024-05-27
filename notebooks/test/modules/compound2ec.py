import requests
import re
from requests.adapters import HTTPAdapter, Retry
import pandas as pd
from Bio.KEGG import REST



# Initializing retry configuration for HTTP requests
retries = Retry(
    total=5,
    backoff_factor=0.25,
    status_forcelist=[500, 502, 503, 504]
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

### KEGG ### 

def compound2ec(id_list: list):
    # Initialize an empty list to store parsed results
    compound2ec_df = []

    # Create a requests session
    session = requests.Session()

    # Iterate over compound IDs
    for compound_id in id_list:
        KEGG_url = f"https://rest.kegg.jp/get/compound:{compound_id}"
        try:
            response = session.get(KEGG_url)
            response.raise_for_status()  # Raise an HTTPError for bad responses

            # Parse the response using Biopython's REST module
            compound_data = REST.kegg_get(compound_id).read()

            # Extract relevant information
            enzymes = None

            for line in compound_data.split("\n"):
                if line.startswith("ENZYME"):
                    enzymes = line[12:].strip()  # Extracting the associated enzyme

            # Append extracted information to results list
            compound2ec_df.append(
                pd.Series({
                    "Entry": compound_id,
                    "Enzyme": enzymes,
                }).to_frame().T
            )
            #print(f'{compound_id} processed')
        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred for {compound_id}: {http_err}")
        except Exception as err:
            print(f"Other error occurred for {compound_id}: {err}")

    # Concatenate the list of parsed results into a DataFrame
    if compound2ec_df:
        compound2ec_df = pd.concat(compound2ec_df, axis=0, ignore_index=True)
    else:
        compound2ec_df = pd.DataFrame(columns=["Entry", "Enzyme"])

    return compound2ec_df

### UniProt ###
