#!/usr/bin/python
import argparse

from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.distributions.OnlineDistribution import OnlineDistribution


class demo_online_encode:
    @staticmethod
    def encode(file, error_correction=nocode, asdna=True):
        epsilon = 0.068
        dist = OnlineDistribution(epsilon)
        number_of_chunks = dist.get_size()
        quality = 7
        dna_rules = FastDNARules()
        if asdna:
            rules = dna_rules
        else:
            rules = None
        encoder = OnlineEncoder(
            file, number_of_chunks, dist, epsilon, quality, error_correction=error_correction, quality_len_format="B",
            insert_header=False, check_block_number_len_format="H", number_of_chunks_len_format="H", rules=rules,
            save_number_of_chunks_in_packet=False)  # , pseudo_decoder=pseudo)
        encoder.set_overhead_limit(1.70)
        encoder.encode_file(split_to_multiple_files=True, save_as_dna=asdna)
        conf = {'error_correction': _error_correction, 'repair_symbols': _repair_symbols, 'quality': quality,
                'epsilon': epsilon, 'find_minimum_mode': False, 'seq_seed': False}
        encoder.save_config_file(conf, section_name="Online_" + file)


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
    _file = args.filename
    _repair_symbols = args.repair_symbols
    _error_correction = get_error_correction_encode(args.error_correction, _repair_symbols)
    print("Zu kodierende Datei: " + str(_file))
    demo = demo_online_encode()
    demo.encode(_file, error_correction=_error_correction)
    # input("Press Enter to continue ...")
