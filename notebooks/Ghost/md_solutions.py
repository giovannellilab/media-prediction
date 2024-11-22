from tqdm import tqdm 
import pandas as pd
import requests
from requests import Session
from requests.adapters import HTTPAdapter, Retry

def _get_session() -> Session:
    retries = Retry(
        total=5,
        backoff_factor=0.25,
        status_forcelist=[500, 502, 503, 504]
    )
    session = Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session

def get_solutions(id_list: list) -> pd.DataFrame:
    session = _get_session()

    base_url = 'https://mediadive.dsmz.de/rest/solution/{}'
    composition_data = []

    for id in tqdm(id_list):
        url = base_url.format(id)
        response = session.get(url)

        if response.status_code == 200:
            data = response.json()

            # Fields
            solution_id = []
            solution_name = []
            component = []
            component_id = []
            component_gl = []
            sub_solutions = []
            sub_component = []
            sub_component_id = []
            steps = []

            # Solution info
            solution_id.append(data.get('data', {}).get('id'))
            solution_name.append(data.get('data', {}).get('name'))
            steps.append(data.get('data', {}).get('steps'))

            recipe = data.get('data', {}).get('recipe', []) 
            for item in recipe:

                ## Solution compounds
                if item.get('compound') is not None and item.get('solution') is None:
                    component.append(item.get('compound'))
                    component_id.append(item.get('compound_id'))
                    component_gl.append(item.get('g_l'))

                ## Sub-solutions/compounds
                elif item.get('solution') is not None and item.get('compound') is None:
                    sub_solutions.append(item.get('solution_id'))

                    for sub_solution in sub_solutions:
                        url = base_url.format(sub_solution)
                        response = session.get(url)

                        if response.status_code == 200:
                            data = response.json()
                            recipe = data.get('data', {}).get('recipe', []) 

                            for item in recipe:
                                sub_component.append(item.get('compound'))
                                sub_component_id.append(item.get('compound_id'))

            # Dataframe preparation
            composition_data.append({
                'solution_id': solution_id,
                'solution_name': solution_name,
                'component': component,
                'component_id': component_id,
                'component_gl': component_gl,
                'sub_sol': sub_solutions,
                'sub_component': sub_component,
                'sub_component_id': sub_component_id,
                'steps': steps
            })
            
        else:
            print(f"Request failed with status code: {response.status_code}")

    composition_df = pd.DataFrame(composition_data)
    return composition_df