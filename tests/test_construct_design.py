from responses import RequestsMock
import pytest
import pandas as pd

from modules import construct_design
from modules.construct_design import TargetData, CODON_TABLE, PrimerDirection
import test_sequences


@pytest.fixture
def mock_web():
    with RequestsMock() as responses:
        yield responses


def test_fetch_target_data(mock_web):
    test_uniprot_id = "fake_uniprot_id"

    fake_returned_data = [{"uniprotSequence": "FAKESEQUENCE", "pdbUrl": "FAKEURL"}]

    expected_output = TargetData(
        uniprot_id="fake_uniprot_id",
        uniprot_sequence="FAKESEQUENCE",
        alphafold_db_url="FAKEURL",
    )

    mock_web.get(
        f"https://alphafold.ebi.ac.uk/api/prediction/{test_uniprot_id}",
        json=fake_returned_data,
        status=200,
    )

    output = construct_design.fetch_target_data(uniprot_id=test_uniprot_id)
    assert (
        mock_web.calls[0].request.url
        == f"https://alphafold.ebi.ac.uk/api/prediction/{test_uniprot_id}"
    )
    assert mock_web.calls[0].request.method == "GET"
    assert output == expected_output


@pytest.mark.parametrize(
    ["input_protein_sequence", "expected_reverse_translation"],
    [
        (
            "TEYKLVVVGAGGVGKSALTIQLIQNH",
            "ACTGAATATAAATTAGTTGTTGTTGGTGCTGGTGGTGTTGGTAAATCTGCTTTAACTATTCAATTAATTCAAAATCAT",
        ),
        ("GLNDIFEAQKIEWHE", "GGTTTAAATGATATTTTTGAAGCTCAAAAAATTGAATGGCATGAA"),
        ("NTAREALPRTSEQ", "AATACTGCTCGTGAAGCTTTACCTCGTACTTCTGAACAA"),
    ],
)
def test_reverse_translate(input_protein_sequence, expected_reverse_translation):
    reverse_translation = construct_design.reverse_translate(
        protein_sequence=input_protein_sequence, table=CODON_TABLE
    )
    assert reverse_translation == expected_reverse_translation


@pytest.mark.parametrize(
    ["input_protein_sequence", "input_template", "input_direction", "expected_primer"],
    [
        (
            "GHFDPAKCRYALGMQDRTIPDSDISASSSWSDSTAARHSRLESSDGDGAWCPAGS",
            test_sequences.Q08345_reverse_translate,
            PrimerDirection.fwd,
            "GGTCATTTTGATCCTGCTAAATGTCG",
        ),
        (
            "GHFDPAKCRYALGMQDRTIPDSDISASSSWSDSTAARHSRLESSDGDGAWCPAGS",
            test_sequences.Q08345_reverse_translate,
            PrimerDirection.rev,
            "AGAACCAGCAGGACACCAAGC",
        ),
        (
            "RHSRLESSDGDGAWCPAGSVFPKEEEYLQVDLQRLHLVALVGTQGRHAGGLGKEFSRSYRLRYSRDGRRWMGWKDRWGQEVISGNEDPEGVVLKDLGPPMVARLVRFYPRADRVMSVC",
            test_sequences.Q08345_reverse_translate,
            PrimerDirection.fwd,
            "CGTCATTCTCGTTTAGAATCTTCTGATG",
        ),
        (
            "RHSRLESSDGDGAWCPAGSVFPKEEEYLQVDLQRLHLVALVGTQGRHAGGLGKEFSRSYRLRYSRDGRRWMGWKDRWGQEVISGNEDPEGVVLKDLGPPMVARLVRFYPRADRVMSVC",
            test_sequences.Q08345_reverse_translate,
            PrimerDirection.rev,
            "ACAAACAGACATAACACGATCAGCAC",
        ),
    ],
)
def test_make_primer(
    input_protein_sequence, input_template, input_direction, expected_primer
):
    primer = construct_design.make_primer(
        protein_sequence=input_protein_sequence,
        template=input_template,
        direction=input_direction,
    )
    assert primer == expected_primer


def test_make_primer_failed():
    expected_error_message = "Could not design primer"

    with pytest.raises(ValueError) as ex:
        construct_design.make_primer(
            protein_sequence=test_sequences.Q08345_small_sequence_fragment,
            template=test_sequences.Q08345_reverse_translate_fragment,
            direction=PrimerDirection.fwd,
        )
    assert str(ex.value) == expected_error_message


example_primer_data = {
    "fwd_primer": ["PRIMER1", "PRIMER2", "PRIMER3"],
    "rev_primer": ["PRIMER4", "PRIMER4", "PRIMER5"],
}


@pytest.mark.parametrize(
    ["input_primer_data", "input_direction", "expected_primer_names"],
    [
        (
            example_primer_data,
            PrimerDirection.fwd,
            {
                "fwd_primer": ["PRIMER1", "PRIMER2", "PRIMER3"],
                "fwd_primer_name": [
                    "fwd_primer_001",
                    "fwd_primer_002",
                    "fwd_primer_003",
                ],
            },
        ),
        (
            example_primer_data,
            PrimerDirection.rev,
            {
                "rev_primer": ["PRIMER4", "PRIMER5"],
                "rev_primer_name": ["rev_primer_001", "rev_primer_002"],
            },
        ),
    ],
)
def test_generate_primer_names(
    input_primer_data, input_direction, expected_primer_names
):
    primer_name_dataframe = construct_design.generate_primer_names(
        input_df=pd.DataFrame(input_primer_data), direction=input_direction
    )
    expected_primer_name_dataframe = pd.DataFrame(expected_primer_names)
    pd.testing.assert_frame_equal(primer_name_dataframe, expected_primer_name_dataframe)


example_bad_primer_data = {
    "not_this_column": ["NOTAPRIMER"],
    "not_this_column_either": ["ALSONOTAPRIMER"],
}


@pytest.mark.parametrize(
    ["input_primer_data", "input_direction", "expected_error_message"],
    [
        (
            example_bad_primer_data,
            PrimerDirection.fwd,
            "fwd_primer not found in input dataframe.",
        ),
        (
            example_bad_primer_data,
            PrimerDirection.rev,
            "rev_primer not found in input dataframe.",
        ),
    ],
)
def test_generate_primer_names_missing_column(
    input_primer_data, input_direction, expected_error_message
):
    with pytest.raises(LookupError) as ex:
        construct_design.generate_primer_names(
            input_df=pd.DataFrame(input_primer_data), direction=input_direction
        )
    assert str(ex.value) == expected_error_message
