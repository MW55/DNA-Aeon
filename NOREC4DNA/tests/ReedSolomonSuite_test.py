import os
import shutil
import pytest
import filecmp
from math import ceil

from norec4dna.helper.bin2Quaternary import string2QUATS
from norec4dna.helper.quaternary2Bin import quats_to_bytes
from norec4dna.ReedSolomonSuite import get_file_size, ReedSolomonEncoder, ReedSolomonDecoder

file = "logo.jpg"
file_dec = "RS_logo.jpg.DECODED"
out_dir = "RS_logo.jpg"


# testing the ReedSolomonSuite with different parameters
@pytest.mark.parametrize("overhead", [0.10])
@pytest.mark.parametrize("chunksize", [100])
@pytest.mark.parametrize("headerchunk", [True])
@pytest.mark.parametrize("flip_bases", [1])
def test_suite(overhead, chunksize, headerchunk, flip_bases):
    file_size = get_file_size(file)
    number_of_chunks = ceil(1.0 * file_size / chunksize) + 1 if headerchunk else 0
    encoder = ReedSolomonEncoder(file, number_of_chunks, 0, overhead)
    encoder.save_packets(True)
    # manipulating some packets
    for fi in os.listdir(out_dir):
        with open(out_dir + "/" + fi, "rb+") as f:
            data = f.read()
            dna_data = "".join(string2QUATS(data))
            dna_data = list(dna_data)
            if dna_data:
                for i in range(16, 16 + flip_bases):
                    if dna_data[i] == 'A':
                        dna_data[i] = 'T'
                    elif dna_data[i] == 'T':
                        dna_data[i] = 'G'
                    elif dna_data[i] == 'G':
                        dna_data[i] = 'C'
                    elif dna_data[i] == 'C':
                        dna_data[i] = 'A'
                dna_data = "".join(dna_data)
                dna_data_temp = b""
                for j in range(0, len(dna_data), 4):
                    try:
                        dna_data_temp += quats_to_bytes(dna_data[j:j + 4])
                    except:
                        pass
                f.seek(0)
                f.write(dna_data_temp)
                f.truncate()
                f.close()
    decoder = ReedSolomonDecoder(out_dir)
    decoder.decodeFolder()
    assert os.path.exists(file_dec)
    assert filecmp.cmp(file_dec, file)
    os.remove(file_dec)
    shutil.rmtree(out_dir)
