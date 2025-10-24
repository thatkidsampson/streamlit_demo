from responses import RequestsMock
import pytest

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
