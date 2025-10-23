from dataclasses import dataclass
import requests

# Base URL for AlphaFoldDB API and headers for requests
ALPHAFOLDDB_BASE_URL = "https://alphafold.ebi.ac.uk/api/prediction"
HEADERS = {
    "Accept": "application/json",
}


@dataclass
class TargetData:
    uniprot_id: str
    uniprot_sequence: str
    alphafold_db_url: str

    @property
    def sequence_length(self) -> int:
        return len(self.uniprot_sequence)


def fetch_target_data(uniprot_id: str) -> TargetData:
    """Fetches target data from AlphaFoldDB for a given UniProt ID."""
    response = requests.get(f"{ALPHAFOLDDB_BASE_URL}/{uniprot_id}", headers=HEADERS)
    if response.status_code != 200:
        raise RuntimeError(
            f"Failed to fetch prediction for UNIPROT ID: {uniprot_id}. The AFDB API returned: {response.status_code} {response.reason}"
        )
    # use the first prediction in the response
    first_prediction = response.json()[0]
    database_info = TargetData(
        uniprot_id=uniprot_id,
        uniprot_sequence=first_prediction["uniprotSequence"],
        alphafold_db_url=first_prediction["pdbUrl"],
    )
    return database_info
