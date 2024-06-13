from tqdm import tqdm
from io import StringIO
import pandas as pd
from modules.utils import _get_session

def _format_cofactors(df: pd.DataFrame) -> pd.DataFrame:
    if "Cofactor" not in df.columns or df["Cofactor"].isnull().all():
        df["CofactorExtracted"] = None
        df["CofactorNormalized"] = None
        df["Fe-S cluster"] = False
        df["CofactorFinal"] = None
        return df

    cofactors_df = df["Cofactor"]\
        .str.extractall(pat=r"(COFACTOR: Name=(.+?); Xref=)").iloc[:, -1]\
        .reset_index(drop=False)\
        .drop("match", axis=1)\
        .set_index("level_0")\
        .rename(columns={1: "CofactorExtracted"})

    cofactors_df = pd.merge(
        left=df,
        right=cofactors_df,
        left_index=True,
        right_index=True,
        how='left'
    )

    # Normalize cofactors
    cofactors_df["CofactorNormalized"] = cofactors_df["CofactorExtracted"]\
        .str.replace(pat="\\(.*\\+\\)", repl="", regex=True)\
        .str.replace(pat=" anion| cation", repl="", regex=True)\
        .str.replace(
            pat="iron-sulfur cluster",
            repl="[Fe-S] cluster",
            regex=True
        )\
        .str.replace(
            pat=".*metal.*",
            repl="metal",
            regex=True
        )\
        .str.replace(
            pat="vanadium|vanadate",
            repl="V",
            regex=True
        )\
        .str.replace(
            pat="^Mo-.*|Momolybdopterin.*|Mobis.*|methylCo",
            repl="Mo",
            regex=True
        )\
        .str.replace(
            pat="^Co-.*|cob.*",
            repl="Co",
            regex=True
        )\
        .str.replace(
            pat="adenosylCo",
            repl="Co",
            regex=False
        )\
        .str.replace(
            pat="methylCo",
            repl="Co",
            regex=False
        )\
        .str.replace(
            pat="^Fe-.*|^Fecoproporphyrin.*|.*heme.*",
            repl="Fe",
            regex=True
        )\
        .str.replace(
            pat="^W-.*",
            repl="W",
            regex=True
        )\
        .str.replace(
            pat="chloride",
            repl="Cl",
            regex=True
        )\
        .str.replace(
            pat="coenzyme F430",
            repl="Ni",
            regex=False
        )\
        .str.replace(
            pat="Ni\\(.*",
            repl="Ni",
            regex=True
        )\
        .str.replace(
            pat=".*copper.*",
            repl="Cu",
            regex=True
        )

    # Extract presence of a Fe-S cluster
    cofactors_df["Fe-S cluster"] = cofactors_df["CofactorNormalized"]\
        .str.contains(
            pat="\\d*Fe-.*\\d*S",
            regex=True,
            na=False
        )

    cofactors_df[
        cofactors_df["Fe-S cluster"] == True
    ][["Fe-S cluster", "CofactorNormalized"]].drop_duplicates()

    # Extract each of the metals in the cluster
    cofactors_df["CofactorFinal"] = cofactors_df["CofactorNormalized"]\
        .str.extract("\\[(.*)\\]", expand=False)\
        .str.replace("\\d+", "", regex=True)\
        .str.split("-")

    # Create a row for each of the metals present in the cluster
    cofactors_df = cofactors_df.explode("CofactorFinal")

    cofactors_df["CofactorFinal"] = cofactors_df["CofactorFinal"].fillna(
        cofactors_df["CofactorNormalized"]
    )

    return cofactors_df

def ec2metals(id_list: list) -> pd.DataFrame:
    fields = [
        "accession",
        "id",
        "gene_names",
        "organism_name",
        "organism_id",
        "lineage",
        "lineage_ids",
        "protein_name",
        "ft_binding",
        "cc_catalytic_activity",
        "cc_cofactor",
        "cc_function",
        "cc_pathway",
        "protein_families",
        "xref_refseq",
        "xref_pdb",
        "xref_kegg",
        "xref_eggnog",
        "xref_bindingdb",
        "xref_chembl",
        "xref_cazy",
        "xref_biocyc",
        "xref_brenda",
        "xref_pathwaycommons",
        "xref_reactome",
        "xref_interpro",
        "xref_pfam"
    ]
    fields = ",".join(fields)

    record_columns = [
        "Query EC",
        "Entry",
        "Entry Name",
        "Gene Names",
        "Organism",
        "Organism (ID)",
        "Taxonomic lineage",
        "Taxonomic lineage (Ids)",
        "Protein names",
        "Binding site",
        "Catalytic activity",
        "Cofactor",
        "Function [CC]",
        "Pathway",
        "Protein families",
        "RefSeq",
        "PDB",
        "KEGG",
        "eggNOG",
        "BindingDB",
        "ChEMBL",
        "CAZy",
        "BioCyc",
        "BRENDA",
        "PathwayCommons",
        "Reactome",
        "InterPro",
        "Pfam"
    ]

    session = _get_session()
    base_url = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=tsv&query=(((keyword%3AKW-0479))+AND+(reviewed%3Atrue)+AND+(ec%3A{}))&fields={}"

    results_list = []

    for ec_number in tqdm(id_list):
        url = base_url.format(ec_number, fields)

        response = session.get(url)
        response.raise_for_status()

        if not len(response.text):
            # Store records not found
            response_df = pd.DataFrame(columns=record_columns, index=range(1))

        else:
            # Read response as a table
            response_df = pd.read_table(StringIO(response.text))

        # Create column for story queried EC (for further plots)
        response_df["Query EC"] = ec_number

        # Sort columns
        response_df = response_df[record_columns].copy()


        # Extract and format cofactors
        response_df = _format_cofactors(response_df)

        results_list.append(response_df)

    return pd.concat(
        results_list,
        axis=0,
        ignore_index=True
    )
