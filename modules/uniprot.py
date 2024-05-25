from tqdm import tqdm

import pandas as pd

from modules.utils import _get_session


def _get_ec_values(record: dict) -> list:
    return "|".join([
        ec_record["value"] for ec_record in
        record\
            .get("proteinDescription", dict())\
            .get("recommendedName", dict())\
            .get("ecNumbers", "")
    ])


def _get_record(record: dict) -> pd.DataFrame:
    record_df = pd.Series({
        "entryType": record["entryType"],
        "primaryAccession": record["primaryAccession"],
        "uniProtkbId": record["uniProtkbId"],
        "taxonId": record["organism"]["taxonId"],
        "fullName": record\
            .get("proteinDescription", dict())\
            .get("recommendedName", dict())\
            .get("fullName", dict())\
            .get("value", None),
        "ecNumbers": _get_ec_values(record)
    })
    record_df = record_df.to_frame().T

    return record_df


def taxon2ec(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = "https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)+AND+(organism_id:{})&fields=id,organism_id,ec,cc_function,cc_pathway"

    results_list = []

    for taxon_id in tqdm(id_list):
        url = base_url.format(taxon_id)

        response = session.get(url)
        response.raise_for_status()

        for record in response.json()["results"]:
            results_list.append(_get_record(record))

    return pd.concat(
        results_list,
        axis=0,
        ignore_index=True
    )
