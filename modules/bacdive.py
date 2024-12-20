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

    taxon_df = taxon_df.reset_index(drop=True)

    # If more than one taxon ID at the species level is present, combine them
    # WARNING: this drops any strain-level ID
    if taxon_df.columns.tolist().count("taxon_id_species") > 1:
        taxon_df = pd.Series({
            "taxon_id_species": ";".join(
                taxon_df["taxon_id_species"].astype(str).values[0]
            )
        }).to_frame().T

    return taxon_df


def taxon2ec(
    id_list: list,
    client: bacdive.client.BacdiveClient
) -> pd.DataFrame:

    client.search(id=id_list)

    results_list = []

    for strain in client.retrieve():
        response_df = pd.json_normalize(strain)

        # -------------------------------------------------------------------- #
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
                "general_ncbi_tax_id": "taxon_id_species",
                "ncbi_tax_id_ncbi_tax_id": "taxon_id_species",
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

        # -------------------------------------------------------------------- #
        # Format taxon IDs

        # Initialize taxon_id_strain column
        response_df["taxon_id_strain"] = None

        # Missing NCBI taxon ID
        if "taxon_id_species" not in response_df.columns:
            response_df["taxon_id_species"] = None

        # Add NCBI taxon ID information
        elif type(response_df["taxon_id_species"].values[0]) == list:
            taxon_df = _format_taxon_id(
                response_df["taxon_id_species"].values[0]
            )
            response_df = pd.concat(
                [
                    taxon_df,
                    response_df.drop(
                        ["taxon_id_species", "taxon_id_strain"],
                        axis=1
                    )
                ],
                axis=1,
                ignore_index=False
            )

        # Rename taxon_id_species to taxon_id
        response_df = response_df.rename(columns={
            "taxon_id_species": "taxon_id"
        })

        # Fix multiple strain IDs for the same record
        if response_df.columns.tolist().count("taxon_id_strain") > 1:
            response_df = response_df.rename(columns={
                "taxon_id_strain": "taxon_id_strain_old"
            })
            response_df["taxon_id_strain"] = ";".join(
                response_df["taxon_id_strain_old"].astype(str).values[0]
            )
            response_df = response_df.drop("taxon_id_strain_old", axis=1)

        # -------------------------------------------------------------------- #
        # Extract EC numbers and store as list

        response_df["ec"] = None

        # If the "ec_numbers" field is present
        if "ec_numbers" in response_df.columns:
            ec_list = pd.DataFrame.from_dict(
                response_df["ec_numbers"].values[0],
                orient="columns"
            )

            # Sometimes the "ec" column is not present
            if "ec" in ec_list.columns:
                response_df["ec"] = [ec_list["ec"].dropna().unique().tolist()]

        # Append to list
        results_list.append(response_df)

    return pd.concat(
        results_list,
        axis=0,
        ignore_index=True
    )
