from tqdm import tqdm

import pandas as pd

import requests

from .utils import _get_session


def get_media() -> pd.DataFrame:

    session = _get_session()
    url = "https://mediadive.dsmz.de/rest/media"

    response = session.get(url)
    response.raise_for_status()

    return pd.DataFrame(response.json()["data"])


def get_composition(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = 'https://mediadive.dsmz.de/rest/medium/{}'

    composition_data = []

    for media_id in tqdm(id_list):
        url = base_url.format(media_id)
        response = session.get(url)

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
                        #component_ids.append(item.get('solution_id')) #don't record solution_id's, only the components
            # Append data for this media to composition_data
            composition_data.append({
                'media_id': media_id,
                'components': components,
                'component_ids': component_ids
            })     
        else:
            print(f"Request failed with status code: {response.status_code}")

    # Convert the list of dictionaries to a DataFrame
    composition_df = pd.DataFrame(composition_data)

    return composition_df


def get_strains(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = 'https://mediadive.dsmz.de/rest/medium-strains/{}'

    strain_data = []

    for media_id in tqdm(id_list):
        url = base_url.format(media_id)
        response = session.get(url)

        if response.status_code == 200:
            data = response.json()
            
            strains = data.get('data', [])
            for strain in strains:
                strain_data.append({
                    'media_id': media_id,
                    'strain_id': strain.get('id'),
                    'species': strain.get('species'),
                    'ccno': strain.get('ccno'),
                    'bacdive_id': strain.get('bacdive_id')
                })

    # Convert the list of dictionaries to a DataFrame
    strain_df = pd.DataFrame(strain_data)

    # Convert BacDive ID to integer
    strain_df["bacdive_id"] = strain_df["bacdive_id"].astype(int)

    return strain_df


def get_compounds(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = 'https://mediadive.dsmz.de/rest/ingredient/{}'

    ingredient_data = []

    for id in tqdm(id_list):
        url = base_url.format(id)
        response = session.get(url)

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


def get_concentrations(id_list: list) -> pd.DataFrame:

    session = _get_session()
    base_url = 'https://mediadive.dsmz.de/rest/medium/{}'

    composition_data = []

    for media_id in tqdm(id_list):
        url = base_url.format(media_id)
        response = session.get(url)

        if response.status_code == 200:
            data = response.json()

            # Initialize lists to store components and ids for each media

            # Names of general solutions
            ingredients = [] 
            ingredient_names = []
            components = []
            component_ids = []
            component_gl = []
            steps = []
            sub_solution = []
            solution_ids = []
            solution_amt = []

            # Traverse the nested structure to extract compounds
            solutions = data.get('data', {}).get('solutions', [])
            for solution in solutions:
                
                ingredients.append(solution.get('id'))
                ingredient_names.append(solution.get('name'))
                recipe = solution.get('recipe', [])
                steps.append(solution.get('steps'))

                for item in recipe:
                    if 'compound' in item:  # Check for 'compound' and 'compound_id'
                        components.append(item.get('compound'))
                        component_ids.append(item.get('compound_id'))
                        component_gl.append(item.get('g_l'))

                    elif 'solution' in item:  # Check for 'solution' and 'solution_id'
                        sub_solution.append(item.get('solution'))
                        solution_ids.append(item.get('solution_id'))
                        solution_amt.append(item.get('amount'))

            # Append data for this media to composition_data
            composition_data.append({
                'media_id': media_id,
                'solutions': ingredients,
                'solution_names': ingredient_names,
                'components': components,
                'component_ids': component_ids,
                'component_gL': component_gl,
                'steps': steps,
                'sub_solutions': sub_solution,
                'solution_ids': solution_ids,
                'solution_ml': solution_amt
            })     
        else:
            print(f"Request failed with status code: {response.status_code}")

    # Convert the list of dictionaries to a DataFrame
    composition_df = pd.DataFrame(composition_data)

    return composition_df
