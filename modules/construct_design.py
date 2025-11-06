from dataclasses import dataclass
import requests
from Bio.Data import CodonTable
from enum import StrEnum
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
from string import ascii_uppercase
from typing import Any

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

# Headers for the input file required by the Merck primer ordering system
MERCK_PRIMER_PLATE_HEADERS = [
    "Plate well",
    "Row",
    "Column",
    "Name",
    "5' Mod",
    "Sequence (5' - 3')",
    "3' Mod",
]


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

    @property
    def template_dna_sequence(self) -> Seq:
        """Returns an representative DNA sequence encoding the UniProt sequence."""
        return Seq(
            reverse_translate(protein_sequence=self.uniprot_sequence, table=CODON_TABLE)
        )


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


def reverse_translate(*, protein_sequence: str, table: Any) -> str:
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
    return dna_sequence


def make_primer(
    protein_sequence: str, template: Seq, direction: PrimerDirection
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


def generate_primer_names(
    *, input_df: pd.DataFrame, direction: PrimerDirection
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


def generate_96_platemap() -> list[str]:
    """Generate a list of 96-well plate references."""
    wells = []
    for i in range(8):
        for c in range(1, 13):
            r = ascii_uppercase[i]
            wells.append(f"{r}{str(c).zfill(2)}")
    return wells


def slice_sequence(sequence: str, start_residue: int, end_residue: int) -> str:
    sliced_sequence = sequence[start_residue:end_residue]
    return sliced_sequence


def generate_primer_dataframe(
    *, construct_dictionary: dict, target_data: TargetData
) -> pd.DataFrame:
    # slice the input sequence to generate the various construct sequences
    df = pd.DataFrame.from_dict(
        construct_dictionary, orient="index", columns=["Start residue", "End residue"]
    )
    df["Sequence"] = df.apply(
        lambda x: slice_sequence(
            target_data.uniprot_sequence, x["Start residue"] - 1, x["End residue"] - 1
        ),
        axis=1,
    )

    # generate the forward and reverse primers for each sequence
    for direction in PrimerDirection:
        df[f"{direction}_primer_annealing"] = df["Sequence"].apply(
            make_primer, template=target_data.template_dna_sequence, direction=direction
        )
        # extend these primers to include the necessary BsaI extensions
        df[f"{direction}_primer"] = (
            BSAI_PRIMER_EXTENSIONS[direction] + df[f"{direction}_primer_annealing"]
        )
        # assign the primers auto-generated names
        primer_names = generate_primer_names(input_df=df, direction=direction)
        # merge the dataframes while preserving the left dataframe's index
        df = pd.merge(
            df.reset_index(), primer_names, how="left", on=f"{direction}_primer"
        ).set_index("index")
        df.index.name = None
    # add 96-well plate well references to the dataframe
    wells_96 = generate_96_platemap()
    df["Plate_well"] = wells_96[: len(df)]
    return df


def generate_384_platemap() -> list[str]:
    """Generate a list of 384-well plate references."""
    wells = []
    for i in range(16):
        for c in range(1, 25):
            r = ascii_uppercase[i]
            wells.append(f"{r}{str(c).zfill(2)}")
    return wells


def make_primer_plate(input_df: pd.DataFrame) -> pd.DataFrame:
    """Create a primer plate dataframe for ordering primers."""
    # generate a 384-well plate map
    wells_384 = generate_384_platemap()
    # create a new dataframe with the 384-well plate well references
    primer_plate = pd.DataFrame(wells_384, columns=["Plate well"])
    primer_plate["Row"] = primer_plate["Plate well"].str[0]
    primer_plate["Column"] = primer_plate["Plate well"].str[1:].astype(int)
    primer_plate.sort_values(by=["Column", "Row"], inplace=True, ascending=True)

    # get the primers from our input dataframe
    primer_sets = []
    for direction in PrimerDirection:
        # get all the fwd primer names and sequences, rename to match the column headers needed by Merck
        primer_set = input_df[
            [f"{direction}_primer_name", f"{direction}_primer"]
        ].copy()
        primer_set.rename(
            {
                f"{direction}_primer_name": "Name",
                f"{direction}_primer": "Sequence (5' - 3')",
            },
            axis=1,
            inplace=True,
        )
        # add the dataframe to the primer sets list
        primer_sets.append(primer_set)
    # combine the two sets of primers into one dataframe of all the fwd + rev primers with sequences
    all_primers = pd.concat(primer_sets)
    # remove the duplicates so it's just a list of unique primers
    unique_primers = all_primers.drop_duplicates(subset="Sequence (5' - 3')").copy()
    unique_primers.reset_index(inplace=True, drop=True)

    # concat to add the unique primers onto the primer_plate dataframe
    primer_plate = pd.concat(
        [primer_plate, unique_primers[["Name", "Sequence (5' - 3')"]]], axis=1
    ).reindex(primer_plate.index)
    # make 3' mod and 5' columns to match the format needed by Merck
    primer_plate["5' Mod"] = ""
    primer_plate["3' Mod"] = ""
    # remove unneeded columns and order the remaining ones to fit the format
    primer_plate = primer_plate[MERCK_PRIMER_PLATE_HEADERS]
    return primer_plate
