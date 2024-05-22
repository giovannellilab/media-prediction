from io import StringIO

import pandas as pd

import requests
from requests.adapters import HTTPAdapter, Retry

from Bio.KEGG.Enzyme import parse


def get_orthologs(text: str) -> pd.DataFrame:

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


def get_media2ec(ec_list: list) -> pd.DataFrame:

    retries = Retry(
        total=5,
        backoff_factor=0.25,
        status_forcelist=[500, 502, 503, 504]
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))


    results_df = []

    for ec_number in ec_list:

        url = f"https://rest.kegg.jp/get/ec:{ec_number}"

        response = session.get(url)
        response.raise_for_status()

        for record in parse(StringIO(response.text)):

            # Get orthologs
            orth_df = get_orthologs(response.text)
            orth_list = ";".join(orth_df["ID"].values)

            results_df.append(
                # TODO: finish!
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


    results_df = pd.concat(
        results_df,
        axis=0,
        ignore_index=True
    )

    return results_df


def get_taxon2ec(id_list: list) -> pd.DataFrame:
    # Initialize results CSV file
    results_path = os.path.join(
        DATA_DIR,
        "kegg",
        "taxon2ec-kegg.csv"
    )
    record_columns = [
        "Gene",
        "EC",
        "KEGG ID"
    ]
    pd.DataFrame(columns=record_columns).to_csv(
        results_path,
        index=False,
        mode="w"
    )

    retries = Retry(
        total=5,
        backoff_factor=0.25,
        status_forcelist=[500, 502, 503, 504]
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    for id_idx, kegg_id in enumerate(id_list):

        # print(f"[+] Retrieving EC numbers for organism: {kegg_id}")

        url = f"https://rest.kegg.jp/link/ec/{kegg_id}"

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

        # Save dataframe to file
        response_df.to_csv(
            results_path,
            index=False,
            header=False,
            mode="a"
        )

        # Delete dataframe to save space
        del response_df

        if id_idx % 10 == 0:
            print(
                f"[+] Processed organism {kegg_id} ({id_idx} / {len(id_list)})"
            )

        time.sleep(0.25)
