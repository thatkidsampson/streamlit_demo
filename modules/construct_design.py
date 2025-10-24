from dataclasses import dataclass
import requests
from Bio.Data import CodonTable

# Base URL for AlphaFoldDB API and headers for requests
ALPHAFOLDDB_BASE_URL = "https://alphafold.ebi.ac.uk/api/prediction"
HEADERS = {
    "Accept": "application/json",
}

# CodonTable for reverse translation
CODON_TABLE = CodonTable.unambiguous_dna_by_id[1]  # Standard table

@dataclass
class TargetData:
    uniprot_id: str
    uniprot_sequence: str
    alphafold_db_url: str

    @property
    def sequence_length(self) -> int:
        return len(self.uniprot_sequence)


def fetch_target_data(*, uniprot_id: str) -> TargetData:
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

def reverse_translate(*, protein_sequence: str, table: CodonTable) -> str:
    """Basic reverse translation of a protein sequence to a DNA sequence."""
    dna_sequence = ""
    for amino_acid in protein_sequence:
        codons = [
            codon for codon, aa in table.forward_table.items() if aa == amino_acid
        ]
        if codons:
            # Just take the first codon for simplicity
            chosen_codon = codons[0]
            dna_sequence += chosen_codon
        else:
            dna_sequence += "NNN"  # Handle unknown amino acids
    return dna_sequence


output = reverse_translate(protein_sequence="NTAREALPRTSEQ", table=CODON_TABLE)
print(output)