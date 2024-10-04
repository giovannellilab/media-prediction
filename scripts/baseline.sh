#!/bin/bash

# Taken from https://unix.stackexchange.com/a/505342
helpFunction()
{
  echo ""
  echo "Usage: $0 -i input_file -o output_file -t num_threads"
  echo -e "\t-i The input sequences file"
  echo -e "\t-o The output CSV file"
  echo -e "\t-t Number of threads to use"
  exit 1 # Exit script after printing help
}

while getopts "i:o:t:" opt
do
  case "$opt" in
    i ) input_file="$OPTARG" ;;
    o ) output_file="$OPTARG" ;;
    t ) num_threads="$OPTARG" ;;
    ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
  esac
done

# Print helpFunction in case parameters are empty
if [ -z "$input_file" ] || [ -z "$output_file" ] || [ -z "$num_threads" ]
then
  echo "Some or all of the parameters are empty";
  helpFunction
fi

# ---------------------------------------------------------------------------- #

sina \
  --in $input_file \
  --out $output_file \
  --db ../data/silva/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb \
  --search-db ../data/silva/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb \
  --threads $num_threads \
  --search \
  --fields full_name,acc,tax_slv,tax_embl_ebi_ena,tax_xref_ncbi \
  --lca-fields tax_slv,tax_embl_ebi_ena,tax_xref_ncbi \
  --search-ignore-super \
  --preserve-order
