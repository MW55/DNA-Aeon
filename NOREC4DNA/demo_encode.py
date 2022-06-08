#!/usr/bin/python
import argparse
from norec4dna import Encoder
from norec4dna.LTEncoder import LTEncoder as LTEncoder
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.rules.DNARules_ErlichZielinski import DNARules_ErlichZielinski
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution

INSERT_HEADER = True
IMPLICIT_MODE = True
NUMBER_OF_CHUNKS_IN_PACKET = False
CHUNK_SIZE = 30


class demo_encode:
    @staticmethod
    def encode(file, error_correction=nocode):
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size=CHUNK_SIZE,
                                                                                 insert_header=INSERT_HEADER)
        print("Number of Chunks=%s" % number_of_chunks)
        dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)
        encoder = LTEncoder(file, number_of_chunks, dist, insert_header=INSERT_HEADER, rules=DNARules_ErlichZielinski(),
                            error_correction=error_correction, number_of_chunks_len_format="H", id_len_format="H",
                            used_packets_len_format="H", save_number_of_chunks_in_packet=NUMBER_OF_CHUNKS_IN_PACKET,
                            implicit_mode=IMPLICIT_MODE)
        encoder.set_overhead_limit(5.00)
        encoder.encode_to_packets()
        encoder.save_packets(split_to_multiple_files=True, save_as_dna=True)
        print("Number of Chunks=%s" % encoder.number_of_chunks)


if __name__ == "__main__":
    # try:
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repairsymbols for ReedSolomon (default=2)")
    args = parser.parse_args()
    filename = args.filename
    _repair_symbols = args.repair_symbols
    _error_correction = get_error_correction_encode(args.error_correction, _repair_symbols)
    print("File to encode: " + str(filename))
    demo = demo_encode()
    demo.encode(filename, error_correction=_error_correction)
    # input("Press Enter to continue ...")
