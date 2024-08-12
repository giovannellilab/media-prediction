from io import StringIO
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from Bio.KEGG.Enzyme import parse
from Bio.KEGG import REST


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


def _get_orthologs(text: str) -> pd.DataFrame:

    # Remove last section
    text = text.split("GENES")[0]

    # Remove previous setion
    text = text.split("ORTHOLOGY")[-1]

    # Replace extra whitespaces
    text = text.strip()

    return pd.read_csv(
        StringIO(text),
        sep="  ",
        names=["ID", "name"]
    )


def media2ec(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = "https://rest.kegg.jp/get/ec:{}"

    results_list = []

    for ec_number in id_list:
        url = base_url.format(ec_number)

        response = session.get(url)
        response.raise_for_status()

        for record in parse(StringIO(response.text)):

            # Get orthologs
            orth_df = _get_orthologs(response.text)
            orth_list = ";".join(orth_df["ID"].values)

            results_list.append(
                pd.Series({
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
                }).to_frame().T
            )

    return pd.concat(
        results_list,
        axis=0,
        ignore_index=True
    )


def taxon2ec(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = "https://rest.kegg.jp/link/ec/{}"

    results_list = []

    for kegg_id in id_list:
        url = base_url.format(kegg_id)

        response = session.get(url)
        response.raise_for_status()

        response_df = pd.read_table(
            StringIO(response.text),
            names=["Gene", "EC"]
        )

        # Remove leading "ec"
        response_df["EC"] = response_df["EC"].str.replace("ec:", "")

        # Add KEGG ID column
        response_df["KEGG ID"] = kegg_id

        results_list.append(response_df)

    return pd.concat(
        results_list,
        axis=0,
        ignore_index=True
    )