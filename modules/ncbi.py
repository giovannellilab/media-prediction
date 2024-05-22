import os

import subprocess

import pandas as pd
from io import StringIO


def get_taxon2ec() -> pd.DataFrame:
    DATA_DIR = "../data/ncbi/"


    annotation_df = pd.read_csv(
        os.path.join(
            DATA_DIR,
            "annotation_metadata.csv"
        )
    )


    # Remove duplicates to avoid redundant calls
    # NOTE: in the future, force for the same accession and taxon in the query
    id_list = annotation_df["Accession"].sample(1000).unique()

    id_list = ",".join(id_list)

    # ------------------------------------------------------------------------ #

    # For example queries: https://www.ncbi.nlm.nih.gov/books/NBK179288/#_chapter6_Examples_

    # Check this answer: https://bioinformatics.stackexchange.com/a/5358
    efetch_cmd = f"efetch -db protein -id {id_list} -format gpc -mode xml"

    # Because multiple EC numbers can be present separated by tab, we need to 
    # change the xtract command. Running the command in the terminal, returns
    # the expanded command: "xtract -insd Protein EC_number" converts into
    # xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -LCUS INSDSeq_locus -SEQ INSDSeq_sequence -group INSDFeature -if INSDFeature_key -equals Protein -clr -pfx "\n" -first "&ACCN,&LCUS" -block INSDQualifier -if INSDQualifier_name -equals EC_number -element INSDQualifier_value -block INSDFeature -unless INSDQualifier_name -equals EC_number -lbl "\-"
    # And now we can add the sep argument (otherwise it is an illegal argument)

    xtract_cmd = "xtract -insd Protein EC_number"


    # Split to fit subprocess format
    efetch_cmd = efetch_cmd.split(" ")
    xtract_cmd = xtract_cmd.split(" ")


    # For the correct use of pipes: https://stackoverflow.com/a/13332300
    # For avoiding memory issues: https://stackoverflow.com/questions/13332268/how-to-use-subprocess-command-with-pipes#comment133789553_64796337
    efetch_res = subprocess.Popen(
        efetch_cmd,
        stdout=subprocess.PIPE
    )
    xtract_res = subprocess.check_output(
        xtract_cmd,
        stdin=efetch_res.stdout
    )
    efetch_res.wait()

    # ------------------------------------------------------------------------ #

    # Decode and split in lines, remove last line
    xtract_res = xtract_res.decode("utf-8").split("\n")[:-1]

    # Set rest of elements as EC numbers separated by comma
    group_ec = lambda x: "\t".join([x[0], ";".join(x[1:])])

    # Fix multiple EC numbers per accession
    xtract_res = [group_ec(line.split("\t")) for line in xtract_res]

    # Results are provided as unicode strings
    results_df = pd.read_table(
        filepath_or_buffer=StringIO("\n".join(xtract_res)),
        names=["Accession", "EC number"]
    )

    return results_df
