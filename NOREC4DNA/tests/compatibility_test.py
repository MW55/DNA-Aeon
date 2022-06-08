import os
import shutil
import pytest
import filecmp
from zipfile import ZipFile

from norec4dna.ErrorCorrection import nocode
from norec4dna import RU10Decoder, LTDecoder, OnlineDecoder
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution

file = "test_files/logo_compatibility.jpg"
file_lst = ["DEC_RU10_logo.jpg", "DEC_LT_logo.jpg", "DEC_ONLINE_logo.jpg"]


@pytest.fixture(scope="session", autouse=True)
def round_around_tests():
    with ZipFile("tests/test_files.zip", 'r') as _zip_file:
        _zip_file.extractall()
    yield
    for f in file_lst:
        if os.path.exists(f):
            os.remove(f)
    if os.path.exists("test_files"):
        shutil.rmtree('test_files')


@pytest.mark.parametrize("args", [("test_files/RU10_logo.jpg", False, 720), ("test_files/RU10_logoH.jpg", True, 721)])
def test_ru10(args):
    x = RU10Decoder(args[0], use_headerchunk=args[1], error_correction=nocode,
                    static_number_of_chunks=args[2])
    x.decode(id_len_format="H", number_of_chunks_len_format="B")
    x.saveDecodedFile(null_is_terminator=False, print_to_output=False)
    if not args[1]:
        out_file = "DEC_RU10_logo.jpg"
    else:
        out_file = file_lst[0]
    remove_nullbytes(out_file)
    assert os.path.exists(out_file) and filecmp.cmp(file, out_file)


@pytest.mark.parametrize("args",
                         [("test_files/ONLINE_logo.jpg", False, 231), ("test_files/ONLINE_logoH.jpg", True, 231)])
def test_online(args):
    decoder = OnlineDecoder(args[0], error_correction=nocode, use_headerchunk=args[1],
                            static_number_of_chunks=args[2])
    decoder.decode(quality_len_format="B", check_block_number_len_format="H",
                   number_of_chunks_len_format="H")
    decoder.saveDecodedFile(null_is_terminator=False, print_to_output=False)
    if not args[1]:
        out_file = "DEC_ONLINE_logo.jpg"
    else:
        out_file = file_lst[2]
    remove_nullbytes(out_file)
    assert os.path.exists(out_file) and filecmp.cmp(file, out_file)


@pytest.mark.parametrize("args", [("test_files/LT_logo.jpg", False, 720), ("test_files/LT_logoH.jpg", True, 721)])
def test_lt(args):
    dist = ErlichZielinskiRobustSolitonDistribution(args[2], seed=2)
    decoder = LTDecoder(args[0], error_correction=nocode, use_headerchunk=args[1],
                        static_number_of_chunks=args[2], implicit_mode=False, dist=dist)
    decoder.decode(number_of_chunks_len_format="H", seed_len_format="I", degree_len_format="H")
    decoder.saveDecodedFile(null_is_terminator=False, print_to_output=False)
    if not args[1]:
        out_file = "DEC_LT_logo.jpg"
    else:
        out_file = file_lst[1]
    remove_nullbytes(out_file)
    assert os.path.exists(out_file) and filecmp.cmp(file, out_file)


def remove_nullbytes(file_path):
    with open(file_path, 'rb+') as f:
        data = f.read()
        data = data.rstrip(b'\x00')
        f.seek(0)
        f.write(data)
        f.truncate()


if __name__ == "__main__":
    with ZipFile("test_files.zip", 'r') as zip_f:
        for name in zip_f.namelist():
            member = zip_f.open(name)
            with open(os.path.basename(name), 'wb') as out_f:
                shutil.copyfileobj(member, out_f)
