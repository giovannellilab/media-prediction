from requests import Session
from requests.adapters import HTTPAdapter, Retry
import os
import re
import pandas as pd
import ast
import glob

def _get_session() -> Session:

    retries = Retry(
        total=5,
        backoff_factor=0.25,
        status_forcelist=[500, 502, 503, 504]
    )
    session = Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    return session

def extract_ec_number(ec_string):
    match = re.search(r'EC:(\d+\.\d+\.\d+\.\d+)', ec_string)
    if match:
        return match.group(1)  # Returns the EC number (e.g., '2.1.1.297')
    return 'N/A'

def parse_gff(gff_file):
    genes = []
    features = {}  # To store features by their IDs

    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            attributes = parts[8]
            attr_dict = {}
            for attribute in attributes.split(';'):
                key_value = attribute.split('=')
                if len(key_value) == 2:
                    key, value = key_value
                    key = key.strip()
                    value = value.strip()
                    if key == 'ec':
                        value = extract_ec_number(value)  # Extract numeric part of EC code
                    attr_dict[key] = value

            feature_info = {
                'seqname': parts[0],
                'source': parts[1],
                'feature': parts[2],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'score': parts[5],
                'strand': parts[6],
                'frame': parts[7],
                'attributes': attr_dict
            }
            
            # Store the feature by its ID
            feature_id = attr_dict.get('ID')
            if feature_id:
                features[feature_id] = feature_info

    # Process features to capture hierarchical relationships
    for feature_id, feature_info in features.items():
        if feature_info['feature'] == 'gene':
            genes.append(feature_info)  # Collect all genes

        # Collect mRNA and CDS features linked to genes
        elif feature_info['feature'] in ['mRNA', 'CDS']:
            parent_id = feature_info['attributes'].get('Parent')
            parent_feature = features.get(parent_id)
            if parent_feature and parent_feature['feature'] == 'gene':
                # Attach mRNA and CDS to their parent gene
                parent_feature.setdefault('children', []).append(feature_info)
                if feature_info['feature'] == 'mRNA':
                    genes.append(feature_info)  # Collect mRNAs

    return genes

def genes_to_dataframe(genes, filename):
    data = {
        'filename': [],
        'seqname': [],
        'source': [],
        'feature': [],
        'ID': [],
        'product': [],
        'ec': [],
        'ko': []
    }
    
    for gene in genes:
        data['filename'].append(filename)
        data['seqname'].append(gene['seqname'])
        data['source'].append(gene['source'])
        data['feature'].append(gene['feature'])
        data['ID'].append(gene['attributes'].get('ID', 'N/A'))
        data['product'].append(gene['attributes'].get('product', 'N/A'))
        data['ec'].append(gene['attributes'].get('ec', 'N/A'))
        data['ko'].append(gene['attributes'].get('ko', 'N/A'))
    
    df = pd.DataFrame(data)
    return df


def process_directory(directory):
    if not os.path.isdir(directory):
        print(f"Directory does not exist: {directory}")
        return None

    # Relaxed pattern to match any files that contain "AssemblySet_DRAM.gff"
    pattern = re.compile(r".*AssemblySet_DRAM\.gff$")
    all_genes = []

    for filename in os.listdir(directory):
        if filename.endswith(".gff") and pattern.match(filename):
            gff_file = os.path.join(directory, filename)
            if not os.path.isfile(gff_file):
                print(f"File does not exist: {gff_file}")
                continue
            print(f"Processing file: {gff_file}")
            genes = parse_gff(gff_file)  # Assuming parse_gff is defined elsewhere
            df_genes = genes_to_dataframe(genes, filename)  # Assuming genes_to_dataframe is defined elsewhere
            all_genes.append(df_genes)

    if not all_genes:
        print("No matching .gff files processed.")
        return None

    final_df = pd.concat(all_genes, ignore_index=True)
    return final_df


def expand_dict_list(df, column):
    # Convert the string representation of the list of dictionaries into actual lists
    df[column] = df[column].apply(ast.literal_eval)
    
    # Explode the column with lists of dictionaries to individual rows
    df_expanded = df.explode(column).reset_index(drop=True)
    
    # Normalize the dictionaries into a flat dataframe
    expanded_rows = pd.json_normalize(df_expanded[column])
    
    # Combine the original dataframe (excluding the original column) with the expanded rows
    result = pd.concat([df_expanded.drop(columns=[column]), expanded_rows], axis=1)
    
    return result