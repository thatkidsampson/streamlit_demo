from responses import RequestsMock
import pytest

from modules import construct_design
from modules.construct_design import TargetData


@pytest.fixture
def mock_web():
    with RequestsMock() as responses:
        yield responses


def test_fetch_target_data(mock_web):
    test_uniprot_id = "fake_uniprot_id"

    fake_returned_data = [{"uniprotSequence": "FAKESEQUENCE",
                           "pdbUrl": "FAKEURL"}]

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

    output = construct_design.fetch_target_data(test_uniprot_id)
    assert (
        mock_web.calls[0].request.url
        == f"https://alphafold.ebi.ac.uk/api/prediction/{test_uniprot_id}"
    )
    assert mock_web.calls[0].request.method == "GET"
    assert output == expected_output
