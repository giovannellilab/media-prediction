#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "No arguments supplied"
    exit 1
fi


while read taxon
do

    echo "[+] Processing taxon ${taxon}"

    # Download protein data for the given taxon ID
    datasets download genome taxon $taxon \
        --annotated \
        --reference \
        --include protein \
        --no-progressbar \
        --filename ${taxon}.zip

    # Unzip into taxon ID folder
    unzip -d $taxon ${taxon}.zip

    # Remove zip file
    rm ${taxon}.zip

    echo

done < $1
