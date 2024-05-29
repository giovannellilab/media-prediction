from io import StringIO
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from Bio.KEGG.Enzyme import parse

def get_orthologs(text: str) -> pd.DataFrame:
    # Remove last section
    text = text.split("GENES")[0]

    # Remove previous section
    text = text.split("ORTHOLOGY")[-1]

    # Replace extra whitespaces
    text = text.strip()

    # Handle empty text case
    if not text:
        return pd.DataFrame(columns=["ID", "name"])

    return pd.read_csv(
        StringIO(text),
        sep="  ",
        names=["ID", "name"],
        engine='python'  # use python engine to handle regex separator
    )

def ec_info(ec_list: list):
    retries = Retry(
        total=5,
        backoff_factor=0.25,
        status_forcelist=[500, 502, 503, 504]
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    results_df = []
    
    for ec_number in ec_list:
        try:
            url = f"https://rest.kegg.jp/get/ec:{ec_number}"
            response = session.get(url)
            response.raise_for_status()

            for record in parse(StringIO(response.text)):
                # Get orthologs
                orth_df = get_orthologs(response.text)
                orth_list = ";".join(orth_df["ID"].values)

                results_df.append(
                    pd.DataFrame([{
                        "Entry": record.entry,
                        "Name": record.name,
                        "Orthologs": orth_list,
                        "Cofactor": record.cofactor,
                        "Pathway": record.pathway,
                        "Reaction": record.reaction,
                        "Product": record.product,
                        "Genes": record.genes,
                        "Structures": record.structures,
                        "Substrate": record.substrate,
                        "DBLinks": record.dblinks
                    }])
                )
        except requests.exceptions.HTTPError as e:
            print(f"HTTP error for EC number {ec_number}: {e}")
        except Exception as e:
            print(f"Error processing EC number {ec_number}: {e}")
            
       

    return results_df



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

