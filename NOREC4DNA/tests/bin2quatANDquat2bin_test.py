from norec4dna.helper import quaternary2Bin, bin2Quaternary
import pytest


@pytest.mark.parametrize("params", [('ACTG'), ('AAAA'), ('CCCC')])
def test_quat2bin_bin2quat(params):
    assert params == bin2Quaternary.byte2QUATS(quaternary2Bin.quats_to_bytes(params))
