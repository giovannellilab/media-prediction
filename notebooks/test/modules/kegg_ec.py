from io import StringIO
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from Bio.KEGG.Enzyme import parse

def get_orthologs(text: str) -> pd.DataFrame:
    # Remove last section
    text = text.split("GENES")[0]

    # Remove previous section
    text = text.split("ORTHOLOGY")[-1]

    # Replace extra whitespaces
    text = text.strip()

    # Handle empty text case
    if not text:
        return pd.DataFrame(columns=["ID", "name"])

    return pd.read_csv(
        StringIO(text),
        sep="  ",
        names=["ID", "name"],
        engine='python'  # use python engine to handle regex separator
    )

def ec_info(ec_list: list):
    retries = Retry(
        total=5,
        backoff_factor=0.25,
        status_forcelist=[500, 502, 503, 504]
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    results_df = []
    
    for ec_number in ec_list:
        try:
            url = f"https://rest.kegg.jp/get/ec:{ec_number}"
            response = session.get(url)
            response.raise_for_status()

            for record in parse(StringIO(response.text)):
                # Get orthologs
                orth_df = get_orthologs(response.text)
                orth_list = ";".join(orth_df["ID"].values)

                results_df.append(
                    pd.DataFrame([{
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
                    }])
                )
        except requests.exceptions.HTTPError as e:
            print(f"HTTP error for EC number {ec_number}: {e}")
        except Exception as e:
            print(f"Error processing EC number {ec_number}: {e}")
            
       

    return results_df

