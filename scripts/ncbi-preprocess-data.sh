#!/bin/bash

# Folder architecture: <TAXON>/ncbi_dataset/data/<ASSEMBLY>/protein.faa

if [[ $# -eq 0 ]] ; then
    echo "No arguments supplied"
    exit 1
fi

# Initialize output file
echo -n "" > all-proteins.faa

for dir in $1/*/
do
    # Get sample name from directory
    taxon=$(basename $dir)

    # See: https://unix.stackexchange.com/a/292254
    # 1) Find all protein.faa inside each folder
    # 2) Add the taxon ID in the header with sed
    # 3) Append to mediadive.faa
    find $dir \
        -name "protein.faa" \
        -exec sed -e "s/>/>${taxon}_/g" {} + >> all-proteins.faa

done
