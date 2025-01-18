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
