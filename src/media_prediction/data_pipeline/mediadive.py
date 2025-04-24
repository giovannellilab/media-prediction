import requests

from tqdm import tqdm

import numpy as np
import pandas as pd

from src.media_prediction.data_pipeline.utils import _get_session


def get_media() -> pd.DataFrame:

    session = _get_session()
    url = "https://mediadive.dsmz.de/rest/media"

    response = session.get(url)
    response.raise_for_status()

    media_df = pd.DataFrame(response.json()["data"])\
        .rename(columns={"id": "media_id"})

    # Convert media_id to string
    media_df["media_id"] = media_df["media_id"].astype(str)

    return media_df


def get_strains(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = "https://mediadive.dsmz.de/rest/medium-strains/{}"

    strain_data = []
    missing_data = []

    for media_id in tqdm(id_list):
        url = base_url.format(media_id)
        response = session.get(url)

        if response.status_code == 200:
            data = response.json()
            
            strains = data.get("data", [])
            for strain in strains:
                strain_data.append({
                    "media_id": media_id,
                    "strain_id": strain.get("id"),
                    "species": strain.get("species"),
                    "ccno": strain.get("ccno"),
                    "bacdive_id": strain.get("bacdive_id")
                })

        else:
            missing_data.append(media_id)

    # Convert the list of dictionaries to a DataFrame
    strain_df = pd.DataFrame(strain_data)

    # Convert media_id to string
    strain_df["media_id"] = strain_df["media_id"].astype(str)

    # Convert BacDive ID to integer
    strain_df["bacdive_id"] = strain_df["bacdive_id"]\
        .fillna(-999)\
        .astype(int)\
        .replace(-999, np.nan)

    print(
        "[WARNING] Failed to retrieve strain data for the following media IDs:",
        missing_data
    )

    return strain_df
