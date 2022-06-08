import os
import shutil
import pytest

from multiplexer.multiplex import Multiplex
from norec4dna.ErrorCorrection import nocode, reed_solomon_encode, crc32

file = '.INFILES/Dorn'
cmp_file = 'tests/cmp_dorn'


# runs before and after every test to cleanup the files
@pytest.fixture(autouse=True)
def run_between_tests():
    if os.path.exists(file):
        os.remove(file)
    shutil.copy(cmp_file, file)


@pytest.mark.parametrize("error_correction", [nocode, reed_solomon_encode, crc32])
@pytest.mark.parametrize("header", [True, False])
@pytest.mark.parametrize("no_channel", [5])
def test_decodable(error_correction, header, no_channel):
    mlt = Multiplex(file, 50, no_channel, 1, error_correction=error_correction, header=header)
    for ind, ch in enumerate(mlt.channel):
        mlt.get_packets_for_channel(200, ind)
    assert mlt.file_potentially_decodable() is True
    for ch in mlt.channel:
        assert ch.is_decodable() is False


@pytest.mark.parametrize("error_correction", [nocode, reed_solomon_encode, crc32])
@pytest.mark.parametrize("header", [True, False])
@pytest.mark.parametrize("no_channel", [3])
@pytest.mark.parametrize("min_shares", [2, 3])
def test_shamirs(error_correction, header, no_channel, min_shares):
    mlt = Multiplex(file, 50, no_channel, 1, error_correction=error_correction, header=header)
    shares = mlt.get_shamir_shares(min_shares)
    data = mlt.combine_shamir_shares(shares)
    assert mlt.no_chunks is data[0]
    assert header is data[1]
    assert error_correction is data[2]


@pytest.mark.parametrize("no_channel", [3])
def test_add_delete_channel(no_channel):
    mlt = Multiplex(file, 50, no_channel, 1)
    assert mlt.remove_channel(0) is False  # Should not be able to delete the only secure channel
    assert mlt.add_channel(False) is True  # Adding public channel
    assert mlt.add_channel(True) is True  # Adding secure channel
    assert mlt.remove_channel(0) is True  # Should be able to delete the first channel now
    assert ("even", 0) not in mlt.used_methods
    mlt.add_channel(True)  # The method "even" should be added with the next channel
    assert ("even", 0) in mlt.used_methods
