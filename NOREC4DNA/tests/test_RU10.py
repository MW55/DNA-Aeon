#!/usr/bin/python
# -*- coding: latin-1 -*-
import os
import shutil
import pytest
import filecmp

from norec4dna import Encoder
from norec4dna import RU10Decoder
from norec4dna import RU10Encoder, RU10BPDecoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.distributions.RaptorDistribution import RaptorDistribution

file = "logo.jpg"
out_dir = "RU10_logo.jpg"
cmp_file = "tests/cmp_logo.jpg"


@pytest.mark.parametrize("as_dna", [False, True])
@pytest.mark.parametrize("decoder_instance", [RU10Decoder])  # :, RU10BPDecoder])
def test_suite(as_dna, decoder_instance):
    dir_path = os.getcwd()
    try:
        os.remove(dir_path + "/" + file)
    except:
        print("Not deleting, File did not exists")
    shutil.copyfile(dir_path + "/" + cmp_file, dir_path + "/" + file)
    print(as_dna)
    chunksize = 200
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = FastDNARules() if as_dna else None
    encoder = RU10Encoder(file, number_of_chunks, dist, pseudo_decoder=pseudo_decoder, rules=rules, id_len_format="H",
                          number_of_chunks_len_format="H", insert_header=True)
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
            pseudo_decoder.is_decoded()
            and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir)
    decoder = decoder_instance(out_dir)
    decoder.decodeFolder(id_len_format="H", number_of_chunks_len_format="H")
    if isinstance(decoder, RU10BPDecoder):
        for pack in encoder.encodedPackets:
            decoder.input_new_packet(pack)
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file)
    decoder.saveDecodedFile(print_to_output=False)
    assert os.path.exists(file) and filecmp.cmp(file, cmp_file)
    shutil.rmtree(out_dir)


if __name__ == "__main__":
    test_suite(False, RU10Decoder)
