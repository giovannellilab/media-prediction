import requests
import pandas as pd

### [MEDIA DSMZ ID TO MEDIA COMPOSITION] ###

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
                    if isinstance(item, dict):  # Check if the item is a dictionary
                        if 'compound' in item:  # Check for 'compound' and 'compound_id'
                            components.append(item.get('compound'))
                            component_ids.append(item.get('compound_id'))
                        elif 'solution' in item:  # Check for 'solution'
                            # Ensure the item is a dictionary
                            solution_data = item.get('solution', {})
                            if isinstance(solution_data, dict):
                                # Extract compounds and their IDs from the nested structure within the solution
                                solution_compounds = solution_data.get('recipe', [])
                                for solution_item in solution_compounds:
                                    if isinstance(solution_item, dict):  # Check if the solution item is a dictionary
                                        if 'compound' in solution_item:
                                            components.append(solution_item.get('compound'))
                                            component_ids.append(solution_item.get('compound_id'))
            # Append data for this media to composition_data
            composition_data.append({
                'media_id': id,
                'components': components,
                'component_ids': component_ids
            })     
        else:
            print(f"Request failed with status code: {response.status_code}")

        #print(f'Retrieved data for {id}')

    # Convert the list of dictionaries to a DataFrame
    composition_df = pd.DataFrame(composition_data)
    return composition_df

### [MEDIA DSMZ ID TO ASSOCIATED TAXA] ###

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
    
        #print(f'Retrieved data for {id}')

    # Convert the list of dictionaries to a DataFrame
    strain_df = pd.DataFrame(strain_data)
    return strain_df

### [MEDIA COMPOSITION TO MEDIA COMPONENTS] ###

def get_compounds(id_list: list):
    base_url = 'https://mediadive.dsmz.de/rest/ingredient/{}'
    ingredient_data = []

    for id in id_list:
        url = base_url.format(id)
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            # Extract data from the 'data' field
            info = data.get('data', {})
            
            # Append data for this component to ingredient_data
            ingredient_data.append({
                'component_id': id,
                'ChEBI': info.get('ChEBI'),
                'KEGG cpd': info.get('KEGG-Compound')
            })
        else:
            print(f"Request for {id} failed with status code: {response.status_code}")
    
        #print(f'Retrieved data for {id}')

    # Convert the list of dictionaries to a DataFrame
    ingr_data = pd.DataFrame(ingredient_data)
    return ingr_data