from dataclasses import dataclass
import requests
from Bio.Data import CodonTable
from Bio import SeqRecord
from enum import StrEnum
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd

# Some parameters for the designed primers:

# the required Tm of the annealing portion of the primers:
PRIMER_TARGET_TM = 60
# BsaI restriction site extensions for the primers:
BSAI_PRIMER_EXTENSIONS = {"fwd": "TATGGTCTCACGAG", "rev": "TATGGTCTCAATGGCTA"}
# When designing primers, we prefer them to end in G or C to reduce mispriming:
PREFERRED_END_NUCLEOTIDES = ["G", "C"]
# minimium primer length:
MIN_PRIMER_LENGTH = 20
# Codon length, in base pairs
CODON_LENGTH = 3

# Base URL for AlphaFoldDB API and headers for requests
ALPHAFOLDDB_BASE_URL = "https://alphafold.ebi.ac.uk/api/prediction"
HEADERS = {
    "Accept": "application/json",
}

# CodonTable for reverse translation
CODON_TABLE = CodonTable.unambiguous_dna_by_id[1]  # Standard table


# StrEnum for primer directions
class PrimerDirection(StrEnum):
    fwd = "fwd"
    rev = "rev"


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


def make_primer(
    protein_sequence: str, template: SeqRecord, direction: PrimerDirection
) -> str:
    """Generates a primer sequence for a given protein sequence and template DNA sequence.
    Args:
        protein_sequence (str): The protein sequence for which to design the primer.
        template (SeqRecord): The template DNA sequence from which to design the primer.
        direction (PrimerDirection): The direction of the primer (forward or reverse).
    Returns:
        str: The designed primer sequence.
    """
    # translate the template DNA sequence into a protein sequence
    translation = template.translate()
    # for a fwd primer, find the start point of the construct protein sequence in the translated template sequence
    if direction == PrimerDirection.fwd:
        loc = translation.find(protein_sequence) * CODON_LENGTH
    # for a rev primer, find the end point of the construct protein sequence in the translated template sequence
    if direction == PrimerDirection.rev:
        loc = (
            translation.find(protein_sequence) + len(protein_sequence)
        ) * CODON_LENGTH
    # start with a primer length of 20bp
    n = MIN_PRIMER_LENGTH
    # start a while loop
    while True:
        # find the sequence n bases from the start or end for a fwd or rev primer
        if direction == PrimerDirection.fwd:
            primer = template[loc : loc + n]
        if direction == PrimerDirection.rev:
            primer = template[loc - n : loc]
            primer = (primer).reverse_complement()
        # get the tm of that primer sequence
        tm = round(mt.Tm_NN(primer, nn_table=mt.DNA_NN2), 2)
        # if the primer has a Tm above the cutoff specified at the top of the notebook and ends in G or C then accept it
        if tm > PRIMER_TARGET_TM and primer[-1] in PREFERRED_END_NUCLEOTIDES:
            return str(primer)
        n += 1
        if n > len(template):
            raise ValueError("Could not design primer")

def generate_primer_names(*,
    input_df: pd.DataFrame, direction: PrimerDirection
) -> pd.DataFrame:
    """Generate auto-names for the primers based on their direction and order."""
    target_column = f"{direction}_primer"
    if target_column not in input_df.columns:
        raise LookupError(f"{target_column} not found in input dataframe.")
    primers_names = input_df[[f"{direction}_primer"]].copy()
    primers_names.drop_duplicates(inplace=True)
    primers_names.reset_index(inplace=True, drop=True)
    primers_names[f"{direction}_primer_name"] = f"{direction}_primer_" + (
        primers_names.index + 1
    ).astype(str).str.zfill(3)
    return primers_names
