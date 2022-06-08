#!/usr/bin/python
import argparse

from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.ErrorCorrection import nocode, get_error_correction_decode
from norec4dna.helper import find_ceil_power_of_four, cluster_and_remove_index, merge_parts, \
    fasta_cluster_and_remove_index

STATIC_NUM_CHUNKS = None  # 149
ID_LEN_FORMAT = "I"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "L"
READ_ALL_BEFORE_DECODER = True


class demo_decode:
    @staticmethod
    def decode(file, error_correction=nocode, null_is_terminator=False, mode_1_bmp=False,
               number_of_chunks=STATIC_NUM_CHUNKS, use_header_chunk=False, id_len_format=ID_LEN_FORMAT,
               number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT, packet_len_format=PACKET_LEN_FORMAT,
               crc_len_format=CRC_LEN_FORMAT, read_all=READ_ALL_BEFORE_DECODER):
        print("Pure Gauss-Mode")
        x = RU10Decoder(file, use_headerchunk=use_header_chunk, error_correction=error_correction,
                        static_number_of_chunks=number_of_chunks)
        x.read_all_before_decode = read_all
        x.decode(id_len_format=id_len_format,
                 number_of_chunks_len_format=number_of_chunks_len_format, packet_len_format=packet_len_format,
                 crc_len_format=crc_len_format)
        x.solve(partial=True)
        if mode_1_bmp:
            return x.mode_1_bmp_decode()
        else:
            return x.saveDecodedFile(null_is_terminator=null_is_terminator, print_to_output=False,
                                     return_file_name=True, partial_decoding=True)


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("filename", metavar="file", type=str, help="the file / folder to Decode")
        parser.add_argument("--error_correction", metavar="error_correction", type=str, required=False,
                            default="nocode", help="Error Correction Method to use; possible values: \
                                                    nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
        parser.add_argument("--repair_symbols", metavar="repair_symbols", type=int, required=False, default=2,
                            help="number of repair_symbols for ReedSolomon (default=2)")
        parser.add_argument("--as_mode_1_bmp", required=False, action="store_true")
        # parser.add_argument("--merge_splits", metavar="merge_splits", required=False, type=int, default=1,
        #                    help="merge from multiple parts")
        parser.add_argument("--number_of_splits", metavar="number_of_splits", required=False, type=int, default=0,
                            help="(optional) number of parts the file has bin split into")
        parser.add_argument("--split_index_position", metavar="split_index_position", required=False, type=str,
                            default="end", help="position of the split index. can be 'start' or 'end")
        parser.add_argument("--split_index_length", metavar="split_index_length", required=False, type=int, default=0,
                            help="number of bases storing the split index")
        parser.add_argument("--last_split_smaller", required=False, action="store_true",
                            help="If set, the number of chunks for the last split will be reduced by 1")
        parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=False, type=int,
                            default=STATIC_NUM_CHUNKS,
                            help="static number of chunks (only set this if not stored in each packet)")
        parser.add_argument("--is_null_terminated", required=False, action="store_true")
        parser.add_argument("--use_header_chunk", required=False, action="store_true")
        args = parser.parse_args()
        _file = args.filename
        _repair_symbols = args.repair_symbols
        _mode_1_bmp = args.as_mode_1_bmp
        _number_of_chunks = args.number_of_chunks
        _last_split_smaller = args.last_split_smaller
        _split_index_position = args.split_index_position
        _split_index_length = args.split_index_length
        _number_of_splits = args.number_of_splits
        _is_null_terminated = args.is_null_terminated
        _use_header_chunk = args.use_header_chunk
        if _number_of_splits != 0:
            _split_index_length = find_ceil_power_of_four(_number_of_splits)
        _last_split_folder = None
        if _split_index_length != 0:
            if _file.lower().endswith("fasta"):
                folders, _last_split_folder = fasta_cluster_and_remove_index(_split_index_position, _split_index_length,
                                                                             _file)
            else:
                folders, _last_split_folder = cluster_and_remove_index(_split_index_position, _split_index_length,
                                                                       _file)
            if _number_of_splits > 0 and _number_of_splits != len(folders):
                print("[WARNING] Number of Splits given by user differs from number of splits found!")
        else:
            folders = [_file]
        _error_correction = get_error_correction_decode(args.error_correction, _repair_symbols)
        decoded_files = []
        for _file in folders:
            print("File / Folder to decode: " + str(_file))
            demo = demo_decode()
            try:
                decoded_files.append(
                    demo.decode(_file, error_correction=_error_correction, null_is_terminator=_is_null_terminated,
                                mode_1_bmp=_mode_1_bmp, number_of_chunks=_number_of_chunks + (
                            -1 if _file == _last_split_folder and _last_split_smaller else 0),
                                use_header_chunk=_use_header_chunk))
            except:
                pass
        if len(folders) > 1:
            merge_parts(decoded_files, remove_tmp_on_success=True)
    except Exception as e:
        raise e
    # input("Press Enter to continue ...")
