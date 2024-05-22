import pandas as pd

import requests


def get_composition(id_list: list):
    base_url = 'https://mediadive.dsmz.de/rest/medium/{}'
    composition_data = []

    for id in id_list:
        url = base_url.format(id)
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            # Initialize lists to store components and ids for each media
            components = []
            component_ids = []
            # Traverse the nested structure to extract compounds
            solutions = data.get('data', {}).get('solutions', [])
            for solution in solutions:
                recipe = solution.get('recipe', [])
                for item in recipe:
                    if 'compound' in item:  # Check for 'compound' and 'compound_id'
                        components.append(item.get('compound'))
                        component_ids.append(item.get('compound_id'))
                    elif 'solution' in item:  # Check for 'solution' and 'solution_id'
                        components.append(item.get('solution'))
                        component_ids.append(item.get('solution_id'))
            # Append data for this media to composition_data
            composition_data.append({
                'media_id': id,
                'components': components,
                'component_ids': component_ids
            })     
        else:
            print(f"Request failed with status code: {response.status_code}")
    
        print(f'Retrieved data for {id}')

    # Convert the list of dictionaries to a DataFrame
    composition_df = pd.DataFrame(composition_data)
    return composition_df


def get_strains(id_list: list):
    base_url = 'https://mediadive.dsmz.de/rest/medium-strains/{}'
    strain_data = []

    for id in id_list:
        url = base_url.format(id)
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            
            strains = data.get('data', [])
            for strain in strains:
                strain_data.append({
                    'media_id': id,
                    'strain_id': strain.get('id'),
                    'species': strain.get('species'),
                    'ccno': strain.get('ccno'),
                    'bacdive_id': strain.get('bacdive_id')
                })     
        else:
            print(f"Request failed with status code: {response.status_code}")
    
        print(f'Retrieved data for {id}')

    # Convert the list of dictionaries to a DataFrame
    strain_df = pd.DataFrame(strain_data)
    return strain_df
