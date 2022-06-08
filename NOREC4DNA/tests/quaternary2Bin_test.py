import pytest

from norec4dna.helper import quaternary2Bin


@pytest.mark.parametrize("params", [('A', 0b00), ('C', 0b01), ('G', 0b10), ('T', 0b11), ('X', None)])
def test_get_quarter_byte(params):
    assert quaternary2Bin.get_quarter_byte(params[0]) == params[1]
