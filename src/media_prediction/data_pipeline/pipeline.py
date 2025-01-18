import os

import requests

import pandas as pd

from src.media_prediction.data_pipeline import mediadive as md


def get_mediadive(data_dir: str) -> pd.DataFrame:

    # Retrieve all available media from MediaDive
    md_media_df = md.get_media()

    # Create list of media IDs
    media_id_list = md_media_df["media_id"].astype(str).unique()

    # Retrieve media-associated strains
    md_strains_df = md.get_strains(media_id_list)

    # Merge media and strains IDs
    data_df = pd.merge(
        left=md_media_df,
        right=md_strains_df,
        on="media_id",
        how="outer",
        indicator="merge_source"
    )
    data_df["merge_source"] = data_df["merge_source"]\
        .cat.rename_categories({
            "right_only": "media_only",
            "left_only": "strains_only"
        })

    data_df.to_csv(
        os.path.join(
            data_dir,
            "mediadive.csv"
        ),
        sep=";",
        index=False
    )

    return data_df
