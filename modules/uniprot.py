import time

import pandas as pd


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


def get_taxon2ec(id_list: list) -> pd.DataFrame:
    # Initialize results CSV file
    results_path = os.path.join(
        DATA_DIR,
        "uniprot",
        "komodo_taxon_to_uniprot_ec.csv"
    )
    record_columns = [
        "entryType",
        "primaryAccession",
        "uniProtkbId",
        "taxonId",
        "fullName",
        "ecNumbers"
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

    for taxon_idx, taxon_id in enumerate(taxon_id_list):
        url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:{taxon_id}&fields=id,organism_id,ec,cc_function,cc_pathway"

        print(f"[+] Retrieving EC numbers for taxon ID: {taxon_id}")

        response = session.get(url)
        response.raise_for_status()

        response_df = pd.DataFrame(columns=record_columns)

        for record_idx, record in enumerate(response.json()["results"]):

            response_df = pd.concat(
                [response_df, _get_record(record)],
                axis=0,
                ignore_index=True
            )

            if record_idx % 10 == 0:
                print(
                    f"[+] Taxon ID {taxon_id} - " + \
                    f"Processed record {record_idx + 1} / " + \
                    f"{len(response.json()['results'])}"
                )

        # Save dataframe to file
        response_df.to_csv(
            results_path,
            index=False,
            header=False,
            mode="a"
        )

        # Delete dataframe to save space
        del response_df

        if taxon_idx % 100 == 0:
            print(f"[+] Processed taxon ID {taxon_idx + 1} / {len(taxon_id_list)}")

        time.sleep(0.25)
