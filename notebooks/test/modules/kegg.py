from io import StringIO
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from Bio.KEGG.Enzyme import parse
from Bio.KEGG import REST



### KEGG Compound Mapping ### 

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
                    "KEGG cpd": compound_id,
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
        compound2ec_df = pd.DataFrame(columns=["KEGG cpd", "Enzyme"])

    return compound2ec_df


### KEGG Orthology ###

def ec2ko(id_list: list):
    # Initialize an empty list to store parsed results
    ec2ko_df = []

    # Create a requests session
    session = requests.Session()

    # Iterate over compound IDs
    for ec in id_list:
        KEGG_url = f"https://rest.kegg.jp/get/ec:{ec}"
        try:
            response = session.get(KEGG_url)
            response.raise_for_status()  # Raise an HTTPError for bad responses

            # Parse the response using Biopython's REST module
            compound_data = REST.kegg_get(ec).read()

            # Extract relevant information
            orthologs = []

            # Split the response into lines and parse
            lines = compound_data.split("\n")
            capture = False
            for line in lines:
                if line.startswith("ORTHOLOGY"):
                    capture = True
                elif capture and not line.startswith(" "):
                    capture = False
                if capture:
                    orthologs.append(line.strip())

            # Combine ortholog lines into a single string with ';' separator
            orthologs = '; '.join(orthologs).replace("ORTHOLOGY", "").strip()

            # Append extracted information to results list
            ec2ko_df.append(
                pd.Series({
                    "ec": ec,
                    "KO": orthologs,
                })
            )
        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred for {ec}: {http_err}")
        except Exception as err:
            print(f"Other error occurred for {ec}: {err}")

    # Convert the list of Series into a DataFrame
    ec2ko_df = pd.DataFrame(ec2ko_df)

    return ec2ko_df

