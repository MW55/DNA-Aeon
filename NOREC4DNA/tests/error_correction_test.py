import pytest
from norec4dna.ErrorCorrection import get_error_correction_decode, get_error_correction_encode


@pytest.mark.parametrize("args", [("nocode", 0), ("reedsolomon", 2), ("reedsolomon", 4),
                                  ("crc", 0)])  # , ("dna_reedsolomon", 2), ("dna_reedsolomon", 4)])
def test_decodability(args):
    e_correction = args[0]
    repair_symbols = int(args[1])
    cmp_str = b"AAAAAACHZACHEFAKLF(24z98"
    assert get_error_correction_decode(e_correction, repair_symbols)(
        get_error_correction_encode(e_correction, repair_symbols)(cmp_str)) == cmp_str
