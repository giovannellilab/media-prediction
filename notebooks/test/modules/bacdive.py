import numpy as np
import pandas as pd

import bacdive


def _format_col(col: str) -> str:
    return "_".join(col.split(".")[-2:])\
        .lower()\
        .replace(" ", "_")\
        .replace("-", "_")


def _format_taxon_id(taxon_values: list) -> pd.DataFrame:
    taxon_df = pd.DataFrame.from_dict(
        taxon_values,
        orient="columns"
    ).T

    taxon_df.columns = taxon_df.loc["Matching level"]
    taxon_df = taxon_df.loc["NCBI tax id"].to_frame().T
    taxon_df.columns = [f"taxon_id_{col}" for col in taxon_df.columns]

    return taxon_df


def taxon2ec(
    id_list: list,
    client: bacdive.client.BacdiveClient
) -> pd.DataFrame:

    client.search(id=id_list)

    results_list = []

    for strain in client.retrieve():
        response_df = pd.json_normalize(strain)

        # Rename columns
        column_mapping = dict(
            zip(
                response_df.columns,
                [_format_col(col) for col in response_df.columns]
            )
        )
        response_df = response_df\
            .rename(columns=column_mapping)\
            .rename(columns={
                "general_bacdive_id": "bacdive_id",
                "general_dsm_number": "dsmz_id",
                "general_ncbi_tax_id": "taxon_id",
                "ncbi_tax_id_ncbi_tax_id": "taxon_id",
                "physiology_and_metabolism_enzymes": "ec_numbers",
                "physiology_and_metabolism_metabolite_utilization": "metabol_uti",
                "metabolite_production_chebi_id": "chebi_id",
                "metabolite_production_metabolite": "metabol_name",
                "metabolite_production_production": "metabol_production",
                "name_and_taxonomic_classification_domain": "domain",
                "name_and_taxonomic_classification_phylum": "phylum",
                "name_and_taxonomic_classification_class": "class",
                "name_and_taxonomic_classification_order": "order",
                "name_and_taxonomic_classification_family": "family",
                "name_and_taxonomic_classification_genus": "genus",
                "name_and_taxonomic_classification_species": "species",
                "name_and_taxonomic_classification_type_strain": "type_strain"
            })

        # Add NCBI taxon ID information
        if type(response_df["taxon_id"].values[0]) == np.ndarray:
            taxon_df = _format_taxon_id(response_df["taxon_id"].values[0])
            response_df = pd.concat(
                [taxon_df.reset_index(drop=True), response_df],
                axis=1,
                ignore_index=False
            )

        # Extract EC numbers and store as list
        if "ec_numbers" in response_df.columns:
            ec_list = pd.DataFrame.from_dict(
                response_df["ec_numbers"].values[0],
                orient="columns"
            )
            response_df["ec"] = [ec_list["ec"].dropna().unique().tolist()]
        else:
            response_df["ec"] = None

        results_list.append(response_df)

    return pd.concat(
        results_list,
        axis=0,
        ignore_index=True
    )