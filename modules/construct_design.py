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

class ExtendedEnum(StrEnum):

    @classmethod
    def list(cls):
        return list(map(lambda c: c.value, cls))

# Headers for the Echo input file
class EchoHeaders(ExtendedEnum):
    source_plate_barcode = "Source Plate Barcode"
    source_plate_name = "Source Plate Name"
    source_plate_type = "Source Plate Type"
    source_well = "Source Well"
    destination_plate_barcode = "Destination Plate Barcode"
    destination_plate_name = "Destination Plate Name"
    destination_well = "Destination Well"
    transfer_volume = "Transfer Volume"
    sample_name = "Sample Name"

# Headers for the input file required by the Merck primer ordering system
class MerckHeaders(ExtendedEnum):
    plate_well = "Plate well"
    row = "Row"
    column = "Column"
    name_ = "Name"
    mod_5 = "5' Mod"
    sequence = "Sequence (5' - 3')"
    mod_3 = "3' Mod"

# StrEnum for primer directions
class PrimerDirection(StrEnum):
    fwd = "fwd"
    rev = "rev"

# set some constants for the Echo input file
# these values will work in an Echo input file but are designed to be easily 
# modifiable to more informative values in a real use case
PRIMER_PLATE_NAME =  "Primer plate"
PRIMER_PLATE_BARCODE = "primer_plate_barcode"
PRIMER_PLATE_TYPE =  "384LDV_AQ_B2" # Beckman low dead volume 384-well plate for aqueous solutions
DESTINATION_PLATE_BARCODE = "PCR_plate_barcode"
DESTINATION_PLATE_NAME = "PCR plate"

# Echo transfer volume in nanolitres
# Here we assume a 5ul reaction volume, 100 µM primer stock concentration, 
# and a required final primer concentration of 250 nM
PCR_REACTION_VOLUME = 5 # µl
PRIMER_CONCENTRATION = 100 # µM
REQUIRED_PRIMER_CONCENTRATION = 250 # nM
ECHO_TRANSFER_VOLUME_NL = (REQUIRED_PRIMER_CONCENTRATION / 1000) / PRIMER_CONCENTRATION * PCR_REACTION_VOLUME * 1000 # nL

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
    elif direction == PrimerDirection.rev:
        loc = (
            translation.find(protein_sequence) + len(protein_sequence)
        ) * CODON_LENGTH
    # If the primer direction is invalid, raise an error
    else:
        raise ValueError("Invalid primer direction specified.")
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


def make_primer_plate(construct_df: pd.DataFrame) -> pd.DataFrame:
    """Create a primer plate dataframe for ordering primers.
    Args:
        construct_df (DataFrame): Details of the primers required to assemble each construct in the format output by generate_primer_dataframe.
    Returns:
        DataFrame: Details of the primer plate to order in the format required by Merck."""
    # generate a 384-well plate map
    wells_384 = generate_384_platemap()
    # create a new dataframe with the 384-well plate well references
    primer_plate = pd.DataFrame(wells_384, columns=[MerckHeaders.plate_well])
    primer_plate[MerckHeaders.row] = primer_plate[MerckHeaders.plate_well].str[0]
    primer_plate[MerckHeaders.column] = primer_plate[MerckHeaders.plate_well].str[1:].astype(int)
    primer_plate.sort_values(by=[MerckHeaders.column, MerckHeaders.row], inplace=True, ascending=True)

    # get the primers from our input dataframe
    primer_sets = []
    for direction in PrimerDirection:
        # get all the fwd primer names and sequences, rename to match the column headers needed by Merck
        primer_set = construct_df[
            [f"{direction}_primer_name", f"{direction}_primer"]
        ].copy()
        primer_set.rename(
            {
                f"{direction}_primer_name": MerckHeaders.name_,
                f"{direction}_primer": MerckHeaders.sequence,
            },
            axis=1,
            inplace=True,
        )
        # add the dataframe to the primer sets list
        primer_sets.append(primer_set)
    # combine the two sets of primers into one dataframe of all the fwd + rev primers with sequences
    all_primers = pd.concat(primer_sets)
    # remove the duplicates so it's just a list of unique primers
    unique_primers = all_primers.drop_duplicates(subset=MerckHeaders.sequence).copy()
    unique_primers.reset_index(inplace=True, drop=True)

    # concat to add the unique primers onto the primer_plate dataframe
    primer_plate = pd.concat(
        [primer_plate, unique_primers[[MerckHeaders.name_, MerckHeaders.sequence]]], axis=1
    ).reindex(primer_plate.index)
    # make 3' mod and 5' columns to match the format needed by Merck
    primer_plate[MerckHeaders.mod_5] = ""
    primer_plate[MerckHeaders.mod_3] = ""
    # remove unneeded columns and order the remaining ones to fit the format
    primer_plate = primer_plate[MerckHeaders.list()]
    return primer_plate


def make_echo_input_file(construct_df: pd.DataFrame, primer_df: pd.DataFrame) -> pd.DataFrame:
    """Create an Echo input file for transferring primers to a PCR plate.
    Args:
        construct_df (DataFrame): Details of the primers required to assemble each construct in the format output by generate_primer_dataframe.
        primer_df (DataFrame): Details of the primer locations in the source place, in the format output by make_primer_plate.
    Returns:
        DataFrame: Picklist for transfer of primers into a PCR plate in the format required by the Echo software."""
    # get all the sequences and required primers
    primer_sets = []
    for direction in PrimerDirection:
        # get the required primer for each sequence
        primer_set = construct_df[["Plate_well", f"{direction}_primer"]].copy()
        # rename the columns to match the Echo input file format
        primer_set.rename(
            {"Plate_well": EchoHeaders.destination_well, f"{direction}_primer": "Primer"},
            axis=1,
            inplace=True,
        )
        # add the dataframe to the primer sets list
        primer_sets.append(primer_set)
    # combine the two sets of primers into one dataframe of all the fwd + rev primers with sequences
    echo_df = pd.concat(primer_sets, ignore_index=True)
    # merge in the primer plate dataframe to get the source plate locations for the primers
    echo_df = pd.merge(
        echo_df,
        primer_df,
        how="left",
        left_on="Primer",
        right_on=MerckHeaders.sequence,
    )

    # set source plate name, source plate barcode, source plate type and transfer volume
    echo_df[EchoHeaders.source_plate_name] = PRIMER_PLATE_NAME
    echo_df[EchoHeaders.source_plate_barcode] = PRIMER_PLATE_BARCODE
    echo_df[EchoHeaders.transfer_volume] = ECHO_TRANSFER_VOLUME_NL
    echo_df[EchoHeaders.source_plate_type] = PRIMER_PLATE_TYPE
    # set the destination plate barcode and name
    echo_df[EchoHeaders.destination_plate_barcode] = DESTINATION_PLATE_BARCODE
    echo_df[EchoHeaders.destination_plate_name] = DESTINATION_PLATE_NAME
    echo_df.rename(
        {MerckHeaders.name_: EchoHeaders.sample_name, MerckHeaders.plate_well: EchoHeaders.source_well}, axis=1, inplace=True
    )
    echo_df = echo_df[EchoHeaders.list()]
    return echo_df
