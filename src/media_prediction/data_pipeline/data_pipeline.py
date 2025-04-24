import os

import requests

import pandas as pd

from src.media_prediction.data_pipeline import mediadive as md


def run_checks(mediadive_dict: dict) -> None:

    media_df = mediadive_dict["media"]
    strains_df = mediadive_dict["strains"]
    media_strains_df = mediadive_dict["media-strains"]

    # ------------------------------------------------------------------------ #

    n_media = media_df["media_id"].nunique()
    n_strains = strains_df["strain_id"].nunique()

    pairs_present = media_strains_df[["media_id", "strain_id"]]\
        .drop_duplicates()\
        .shape[0]

    pairs_data = media_strains_df[
        media_strains_df["merge_source"] == "both"
    ][["media_id", "strain_id"]].drop_duplicates().shape[0]

    pairs_media = media_strains_df[
        media_strains_df["merge_source"] == "media_only"
    ][["media_id", "strain_id"]].drop_duplicates().shape[0]

    pairs_strain = media_strains_df[
        media_strains_df["merge_source"] == "strain_only"
    ][["media_id", "strain_id"]].drop_duplicates().shape[0]

    # ------------------------------------------------------------------------ #

    print("[+] Number of unique media IDs: ", n_media)
    print("[+] Number of unique strain IDs:", n_strains)
    print(
        "[+] Number of unique media-strain pairs present:\t",
        pairs_present
    )
    print(
        "[+] Number of unique media-strain pairs with media data:",
        pairs_data
    )
    print(
        "[+] Number of unique media-strain pairs with only media data: ",
        pairs_media
    )
    print(
        "[+] Number of unique media-strain pairs with only strain data:",
        pairs_strain
    )

    # ------------------------------------------------------------------------ #

    assert n_media == 3323, "[ERROR] Inconsistent number of media!"
    assert n_strains == 46702, "[ERROR] Inconsistent number of strains!"
    assert pairs_present == 71165, "[ERROR] Inconsistent number of pairs!"
    assert pairs_data == 70534, "[ERROR] Inconsistent number of pairs w/ data!"
    assert pairs_media == 631, "[ERROR] Inconsistent number of unpaired media!"
    assert pairs_strain == 0, "[ERROR] Inconsistent number of unpaired strains!"

    return None


def get_mediadive(data_dir: str) -> dict:

    # Retrieve all available media from MediaDive
    print("[+] Retrieving media data...")
    media_df = md.get_media()

    # Create list of media IDs
    media_id_list = media_df["media_id"].unique()

    # Retrieve media-associated strains
    print("[+] Retrieving strains data...")
    strains_df = md.get_strains(media_id_list)

    # Retrieve media-associated components
    print("[+] Retrieving ingredients data...")
    ingredients_df = md.get_ingredients()
    ingredients_df = md.get_ingredients_metadata(ingredients_df["id"].unique())

    # ------------------------------------------------------------------------ #
    # Merge data from different sources to run the checks

    # Merge media and strains IDs
    medium_strains_df = pd.merge(
        left=media_df,
        right=strains_df,
        on="media_id",
        how="outer",
        indicator="merge_source"
    )
    medium_strains_df["merge_source"] = medium_strains_df["merge_source"]\
        .cat.rename_categories({
            "left_only": "media_only",
            "right_only": "strains_only"
        })

    mediadive_dict = {
        "media": media_df,
        "strains": strains_df,
        "ingredients": ingredients_df,
        "media-strains": medium_strains_df
    }

    run_checks(mediadive_dict)

    # ------------------------------------------------------------------------ #
    # Save MediaDive data

    media_df.to_csv(
        os.path.join(
            data_dir,
            "mediadive-media.csv"
        ),
        sep=";",
        index=False
    )
    strains_df.to_csv(
        os.path.join(
            data_dir,
            "mediadive-medium-strains.csv"
        ),
        sep=";",
        index=False
    )
    ingredients_df.to_csv(
        os.path.join(
            data_dir,
            "mediadive-ingredients.csv"
        ),
        sep=";",
        index=False
    )

    return mediadive_dict
