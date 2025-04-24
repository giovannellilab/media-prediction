import os

import requests

import pandas as pd

from src.media_prediction.data_pipeline import mediadive as md


def get_mediadive(data_dir: str) -> dict:

    # Retrieve all available media from MediaDive
    media_df = md.get_media()

    # Create list of media IDs
    media_id_list = media_df["media_id"].astype(str).unique()

    # Retrieve media-associated strains
    strains_df = md.get_strains(media_id_list)

    # ------------------------------------------------------------------------ #
    # Merge data from different sources

    # Merge media and strains IDs
    media_strains_df = pd.merge(
        left=media_df,
        right=strains_df,
        on="media_id",
        how="outer",
        indicator="merge_source"
    )
    media_strains_df["merge_source"] = media_strains_df["merge_source"]\
        .cat.rename_categories({
            "left_only": "media_only",
            "right_only": "strains_only"
        })

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
            "mediadive-strains.csv"
        ),
        sep=";",
        index=False
    )
    media_strains_df.to_csv(
        os.path.join(
            data_dir,
            "mediadive-media-strains.csv"
        ),
        sep=";",
        index=False
    )

    return {
        "media": media_df,
        "strains": strains_df,
        "media-strains": media_strains_df
    }
