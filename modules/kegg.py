from io import StringIO

import pandas as pd

from Bio.KEGG.Enzyme import parse

from modules.utils import _get_session


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
