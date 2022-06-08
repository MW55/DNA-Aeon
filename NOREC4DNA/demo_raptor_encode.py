#!/usr/bin/python
import os
import argparse

from norec4dna.Encoder import Encoder
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper import split_file, number_to_base_str, find_ceil_power_of_four, merge_folder_content

ID_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"
DEFAULT_CHUNK_SIZE = 40


class demo_raptor_encode:
    @staticmethod
    def encode(file, asdna=True, chunk_size=DEFAULT_CHUNK_SIZE, error_correction=nocode, insert_header=False,
               save_number_of_chunks_in_packet=False, mode_1_bmp=False, prepend="", append="", upper_bound=0.5,
               save_as_fasta=True, save_as_zip=True, overhead=0.40):
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
        dist = RaptorDistribution(number_of_chunks)
        if asdna:
            rules = FastDNARules()
        else:
            rules = None
        x = RU10Encoder(file, number_of_chunks, dist, insert_header=insert_header, pseudo_decoder=None,
                        chunk_size=0, rules=rules, error_correction=error_correction,
                        packet_len_format=PACKET_LEN_FORMAT,
                        crc_len_format=CRC_LEN_FORMAT, number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
                        id_len_format=ID_LEN_FORMAT, save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
                        mode_1_bmp=mode_1_bmp, prepend=prepend, append=append, drop_upper_bound=upper_bound)
        x.set_overhead_limit(overhead)
        x.encode_to_packets()
        if save_as_fasta and asdna:
            x.save_packets_fasta(file_ending="_RU10", seed_is_filename=True)
        elif save_as_zip:
            x.save_packets_zip(save_as_dna=asdna, file_ending="_RU10", seed_is_filename=True)
        else:
            x.save_packets(True, save_as_dna=asdna, seed_is_filename=True, clear_output=True)

        return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--as_dna", help="convert packets to dna and use dna rules", action="store_true",
                        required=False)
    parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=DEFAULT_CHUNK_SIZE,
                        help="size of chunks to split the file into")
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", action="store_true", required=False, default=False)
    parser.add_argument("--save_number_of_chunks", metavar="save_number_of_chunks", required=False, type=bool,
                        default=False)
    parser.add_argument("--save_as_fasta", action="store_true", required=False)
    parser.add_argument("--save_as_zip", action="store_true", required=False)
    parser.add_argument("--as_mode_1_bmp", action="store_true",
                        help="convert to a header-less B/W BMP format. (use only for image/bmp input)", required=False)
    parser.add_argument("--split_input", metavar="split_input", required=False, type=int, default=1,
                        help="number of subcodings to split input file into")
    parser.add_argument("--drop_upper_bound", metavar="drop_upper_bound", required=False, type=float, default=0.5,
                        help="upper bound for calculated error probability of packet before dropping")
    parser.add_argument("--overhead", metavar="overhead", required=False, type=float, default=0.40, help="desired overhead of packets")

    args = parser.parse_args()
    _file = args.filename
    _as_dna = args.as_dna
    _chunk_size = args.chunk_size
    _no_repair_symbols = args.repair_symbols
    _number_of_splits = args.split_input
    _insert_header = args.insert_header
    _save_number_of_chunks = args.save_number_of_chunks
    _mode_1_bmp = args.as_mode_1_bmp
    _upper_bound = args.drop_upper_bound
    _error_correction = get_error_correction_encode(args.error_correction, _no_repair_symbols)
    _save_as_fasta = args.save_as_fasta
    _save_as_zip = args.save_as_zip
    _overhead = args.overhead
    if _number_of_splits > 1:
        input_files = split_file(_file, _number_of_splits)
        power_of_four = find_ceil_power_of_four(len(input_files))
        print("Spltting input into {} sub files. We need to prepend/append {} base(s)".format(len(input_files),
                                                                                              power_of_four))
        prepend_matching = {input_files[i]: number_to_base_str(i, power_of_four) for i in range(len(input_files))}
        print("Matching: {}".format(prepend_matching))
    else:
        input_files = [_file]
        prepend_matching = {_file: ""}
    for _file in input_files:
        print("File to encode: " + str(_file))
        demo = demo_raptor_encode()
        encoder_instance = demo.encode(_file, _as_dna, chunk_size=_chunk_size, error_correction=_error_correction,
                                       mode_1_bmp=_mode_1_bmp, insert_header=_insert_header,
                                       append=prepend_matching[_file],
                                       save_number_of_chunks_in_packet=_save_number_of_chunks, upper_bound=_upper_bound,
                                       save_as_fasta=_save_as_fasta, save_as_zip=_save_as_zip, overhead=_overhead)
        conf = {'error_correction': args.error_correction, 'repair_symbols': _no_repair_symbols, 'asdna': _as_dna,
                'number_of_splits': _number_of_splits}
        config_filename = encoder_instance.save_config_file(conf)
        print("Saved config file: %s" % config_filename)

    if len(input_files) > 1:
        merge_folder_content(os.path.dirname(os.path.realpath(_file)), _file + "combined_split_output",
                             append_folder_name=True, clear_dest_folder=True)
    # input("Press Enter to continue ...")
