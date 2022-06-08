from norec4dna.helper import bin2Quaternary
import pytest


@pytest.mark.parametrize("params", [(0, 0, 'A'), (0, 1, 'C'), (1, 0, 'G'), (1, 1, 'T')])
def test_get_quat(params):
    assert bin2Quaternary.getQUAT(params[0], params[1]) == params[2]


@pytest.mark.parametrize("params", [(0b00011011, 'ACGT'), (0b11111111, 'TTTT'), (0b10101010, 'GGGG')])
def test_byte2_quats(params):
    assert bin2Quaternary.byte2QUATS(params[0]) == params[1]
