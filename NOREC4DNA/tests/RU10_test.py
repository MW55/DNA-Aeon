import os
import shutil
import pytest
import filecmp

from norec4dna import Encoder
from norec4dna import RU10Encoder
from norec4dna.rules.DNARules import DNARules
from norec4dna.rules.DNARules2 import DNARules2
from norec4dna import RU10Decoder, RU10BPDecoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode, crc32, reed_solomon_decode, reed_solomon_encode, crc32_decode

file = "logo.jpg"
file2 = "Dorn"
out_dir = "RU10_logo.jpg"
out_dir2 = "RU10_Dorn"
cmp_file = "tests/cmp_logo.jpg"
cmp_file2 = "tests/cmp_dorn"


# runs before and after every test to cleanup the files
@pytest.fixture(autouse=True)
def run_between_tests():
    if os.path.exists(file):
        os.remove(file)
    shutil.copy(cmp_file, file)
    if os.path.exists(file2):
        os.remove(file2)
    shutil.copy(cmp_file2, file2)


# testing as_dna and headerchunk with both dna_rules
@pytest.mark.parametrize("as_dna", [False, True])
@pytest.mark.parametrize("chunk_size", [100])
@pytest.mark.parametrize("dna_rules", [DNARules2()])
@pytest.mark.parametrize("error_correction", [nocode])
@pytest.mark.parametrize("headerchunk", [True, False])
def test_suite(as_dna, chunk_size, dna_rules, error_correction, headerchunk):
    chunksize = chunk_size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    decoder_instance = RU10Decoder
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks + (1 if headerchunk else 0))
    rules = dna_rules if as_dna else None
    encoder = RU10Encoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules, error_correction=error_correction,
        insert_header=headerchunk
    )
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir)
    decoder = decoder_instance(out_dir, error_correction=error_correction, use_headerchunk=headerchunk)
    decoder.decode()
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file)
    decoder.saveDecodedFile(print_to_output=False)
    if headerchunk:
        assert os.path.exists(file) and filecmp.cmp(file, cmp_file)
    else:
        assert os.path.exists("DEC_RU10_logo.jpg") and filecmp.cmp("DEC_RU10_logo.jpg", cmp_file)
    shutil.rmtree(out_dir)


# testing error_correction methods with both dna_rules
@pytest.mark.parametrize("as_dna", [True])
@pytest.mark.parametrize("chunk_size", [100])
@pytest.mark.parametrize("dna_rules", [DNARules(), DNARules2(), FastDNARules()])
@pytest.mark.parametrize("error_correction_pair",
                         [(nocode, nocode), (crc32, crc32_decode), (reed_solomon_encode, reed_solomon_decode)])
def test_suite2(as_dna, chunk_size, dna_rules, error_correction_pair):
    chunksize = chunk_size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    decoder_instance = RU10Decoder
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = dna_rules if as_dna else None
    encoder = RU10Encoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules,
        error_correction=error_correction_pair[0],
    )
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir)
    decoder = decoder_instance(out_dir, error_correction=error_correction_pair[1])
    decoder.decode()
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file)
    decoder.saveDecodedFile(print_to_output=False)
    assert os.path.exists(file) and filecmp.cmp(file, cmp_file)
    shutil.rmtree(out_dir)


# testing different id_len_formats and number_of_chunks_len_formats
# the Dorn txt file is used for the tests
@pytest.mark.parametrize("as_dna", [True])
@pytest.mark.parametrize("chunk_size", [100])
@pytest.mark.parametrize("dna_rules", [None])
@pytest.mark.parametrize("id_len_form", ["H", "I", "B"])
@pytest.mark.parametrize("number_of_chunks_len_form", ["H", "I", "B"])
def test_suite3(as_dna, chunk_size, dna_rules, id_len_form, number_of_chunks_len_form):
    chunksize = chunk_size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file2, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    decoder_instance = RU10Decoder
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = dna_rules if as_dna else None
    encoder = RU10Encoder(
        file2, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules, id_len_format=id_len_form,
        number_of_chunks_len_format=number_of_chunks_len_form
    )
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir2)
    decoder = decoder_instance(out_dir2)
    decoder.decode(id_len_format=id_len_form, number_of_chunks_len_format=number_of_chunks_len_form)
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file2)
    decoder.saveDecodedFile(print_to_output=False)
    assert os.path.exists(file2) and filecmp.cmp(file2, cmp_file2)
    shutil.rmtree(out_dir2)


# testing for corrupt packets with crc32 and reedSolomon and failed decoding
@pytest.mark.parametrize("as_dna", [True])
@pytest.mark.parametrize("chunk_size", [100])
@pytest.mark.parametrize("dna_rules", [None])
@pytest.mark.parametrize("error_correction", [(crc32, crc32_decode), (reed_solomon_encode, reed_solomon_decode)])
def test_suite4(as_dna, chunk_size, dna_rules, error_correction):
    chunksize = chunk_size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    decoder_instance = RU10Decoder
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = dna_rules if as_dna else None
    encoder = RU10Encoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules, error_correction=error_correction[0])
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir)
    # do not delete all packets (and break the last one).
    # that way the GEPP inside the decoder will get initialized and we might not end in a race-condition for
    # decoder.decode() sometimes raising an Exception..
    for i in range(2, number_of_chunks):
        tmp_path = "RU10_" + file + "/" + str(i) + ".RU10_DNA"
        os.remove(tmp_path)
    with open("RU10_" + file + "/0.RU10_DNA", 'rb+') as tmp_file:
        # TODO we should flip bits in the middle rather than deleting 4 bytes at the end (we store crc32 / reedsolomon at the end)
        tmp_file.seek(-4, os.SEEK_END)
        tmp_file.truncate()
    decoder = decoder_instance(out_dir, error_correction=error_correction[1])
    decoder.decode()
    assert decoder.corrupt == 1
    assert not decoder.is_decoded()
    os.remove(file)
    with pytest.raises(AssertionError):
        decoder.saveDecodedFile(print_to_output=False, partial_decoding=False)
    assert not (os.path.exists(file) and filecmp.cmp(file, cmp_file))
    shutil.rmtree(out_dir)


# testing the null_is_terminator option for a txt file without headerchunk and print_to_output
@pytest.mark.parametrize("as_dna", [True])
@pytest.mark.parametrize("chunk_size", [100])
@pytest.mark.parametrize("dna_rules", [None])
@pytest.mark.parametrize("error_correction", [nocode])
@pytest.mark.parametrize("headerchunk", [False])
@pytest.mark.parametrize("decoder_instance", [RU10Decoder])  # , RU10BPDecoder
def test_suite5(as_dna, chunk_size, dna_rules, error_correction, headerchunk, decoder_instance):
    chunksize = chunk_size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file2, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = dna_rules if as_dna else None
    encoder = RU10Encoder(
        file2, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules, error_correction=error_correction,
        insert_header=headerchunk)
    encoder.encode_to_packets()
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir2)
    decoder = decoder_instance(out_dir2, use_headerchunk=headerchunk, error_correction=error_correction)
    decoder.decode()
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file2)
    decoder.saveDecodedFile(print_to_output=True, null_is_terminator=True)
    assert os.path.exists('DEC_RU10_' + file2) and filecmp.cmp('DEC_RU10_' + file2, cmp_file2)
    shutil.rmtree(out_dir2)


""""# testing the null_is_terminator option for a txt file without headerchunk and print_to_output
@pytest.mark.parametrize("as_dna", [True])
@pytest.mark.parametrize("chunk_size", [100])
@pytest.mark.parametrize("dna_rules", [None])
def test_suite5(as_dna, chunk_size, dna_rules):
    chunksize = chunk_size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file2, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    decoder_instance = RU10Decoder
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = dna_rules if as_dna else None
    encoder = RU10Encoder(
        file2, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules)
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=False, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(file2 + ".RU10_DNA")
    decoder = decoder_instance(file2 + ".RU10_DNA")
    decoder.decode()
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file2)
    decoder.saveDecodedFile(print_to_output=True, null_is_terminator=True)
    assert os.path.exists(file2) and filecmp.cmp(file2, cmp_file2)
    os.remove(file2 + ".RU10_DNA")"""
