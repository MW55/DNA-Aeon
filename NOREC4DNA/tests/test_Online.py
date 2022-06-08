#!/usr/bin/python
# -*- coding: latin-1 -*-
import os
import shutil
import filecmp
import pytest

from norec4dna.Encoder import Encoder
from norec4dna.OnlineDecoder import OnlineDecoder
from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.OnlineBPDecoder import OnlineBPDecoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.distributions.OnlineDistribution import OnlineDistribution

file = "logo.jpg"
out_dir = "ONLINE_logo.jpg"
cmp_file = "tests/cmp_logo.jpg"


# WARNING: THIS TEST MIGHT FAIL FOR OTHER INPUT-FILES OR CHUNKSIZES IF use_header IS False (filecmp detects diff at end)
@pytest.mark.parametrize("as_dna", [False, True])
@pytest.mark.parametrize("decoder_instance", [OnlineDecoder, OnlineBPDecoder])
@pytest.mark.parametrize("use_header", [True, False])
def test_suite(as_dna, decoder_instance, use_header):
    dir_path = os.getcwd()
    try:
        os.remove(dir_path + "/" + file)
    except:
        print("Not deleting, File did not exists")
    shutil.copyfile(dir_path + "/" + cmp_file, dir_path + "/" + file)
    chunksize = 200
    epsilon = 0.07
    quality = 5
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = OnlineDistribution(epsilon)
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = FastDNARules() if as_dna else None
    encoder = OnlineEncoder(file, number_of_chunks, dist, epsilon, quality, number_of_chunks_len_format="H",
                            check_block_number_len_format="H", quality_len_format="B", pseudo_decoder=pseudo_decoder,
                            rules=rules, insert_header=use_header)
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (pseudo_decoder.is_decoded() and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks)
    assert os.path.exists(out_dir)
    decoder = decoder_instance(out_dir, use_headerchunk=use_header)
    decoder.decodeFolder(number_of_chunks_len_format="H", check_block_number_len_format="H", quality_len_format="B")
    if decoder_instance == OnlineBPDecoder:
        decoder.decodeFolder(number_of_chunks_len_format="H", check_block_number_len_format="H", quality_len_format="B")
    decoder.solve()
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    if not use_header:
        out_file = "DEC_ONLINE_" + file
    else:
        out_file = file
    try:
        os.remove(out_file)
    except:
        print("Not deleting, File did not exists")
    decoder.saveDecodedFile(print_to_output=False)
    assert os.path.exists(out_file) and filecmp.cmp(out_file, cmp_file)
    if decoder_instance == OnlineBPDecoder:
        # since ApproxDecoder defines an upper bound Gauss-Decoder MUST be able to decode!
        decoder = OnlineDecoder(out_dir, use_headerchunk=use_header)
        decoder.decodeFolder(number_of_chunks_len_format="H",
                             check_block_number_len_format="H",
                             quality_len_format="B")
        assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
        os.remove(out_file)
        decoder.saveDecodedFile(print_to_output=False)
        assert os.path.exists(out_file) and filecmp.cmp(file, cmp_file)
        shutil.rmtree(out_dir)
