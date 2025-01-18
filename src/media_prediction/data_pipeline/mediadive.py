from tqdm import tqdm

import pandas as pd

import requests

from src.media_prediction.data_pipeline.utils import _get_session


def get_media() -> pd.DataFrame:

    session = _get_session()
    url = "https://mediadive.dsmz.de/rest/media"

    response = session.get(url)
    response.raise_for_status()

    return pd.DataFrame(response.json()["data"])\
        .rename(columns={"id": "media_id"})


def get_strains(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = "https://mediadive.dsmz.de/rest/medium-strains/{}"

    strain_data = []

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

    # Convert the list of dictionaries to a DataFrame
    strain_df = pd.DataFrame(strain_data)

    # Convert BacDive ID to integer
    strain_df["bacdive_id"] = strain_df["bacdive_id"].astype(int)

    return strain_df
