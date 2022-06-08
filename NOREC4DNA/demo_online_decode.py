#!/usr/bin/python
import argparse

from norec4dna.OnlineDecoder import OnlineDecoder
from norec4dna.ErrorCorrection import nocode, get_error_correction_decode

STATIC_NUM_CHUNKS = 231
NULL_IS_TERMINATOR = True
PRINT_TO_OUTPUT = True
ID_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "H"
PACKET_LEN_FORMAT = "H"
CRC_LEN_FORMAT = "L"


class demo_decode:
    @staticmethod
    def decode(file, error_correction=nocode, null_is_terminator=NULL_IS_TERMINATOR, mode_1_bmp=False,
               number_of_chunks=STATIC_NUM_CHUNKS, use_header_chunk=False, id_len_format=ID_LEN_FORMAT,
               number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT, packet_len_format=PACKET_LEN_FORMAT,
               crc_len_format=CRC_LEN_FORMAT):
        def _internal(decoder):
            decoder.decode(quality_len_format="B", check_block_number_len_format=id_len_format,
                           number_of_chunks_len_format=number_of_chunks_len_format, crc_len_format=crc_len_format)

        x = OnlineDecoder(file, error_correction=error_correction, use_headerchunk=use_header_chunk,
                          static_number_of_chunks=number_of_chunks)
        print("[1/2] Approximation Decode")
        _internal(x)
        x.saveDecodedFile(null_is_terminator=null_is_terminator, print_to_output=PRINT_TO_OUTPUT)


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "filename", metavar="file", type=str, help="the file / folder to Decode"
        )
        parser.add_argument("--error_correction", metavar="error_correction", type=str, required=False,
                            default="nocode",
                            help="Error Correction Method to use; possible values: \
                                nocode, crc, reedsolomon (default=nocode)")
        parser.add_argument("--repair_symbols", metavar="repair_symbols", type=int, required=False, default=2,
                            help="number of repairsymbols for ReedSolomon (default=2)")
        args = parser.parse_args()
        _file = args.filename
        _repair_symbols = args.repair_symbols
        _error_correction = get_error_correction_decode(args.error_correction, _repair_symbols)
        print("Zu dekodierende Datei / Ordner: " + str(_file))
        demo = demo_decode()
        demo.decode(_file, error_correction=_error_correction)
        # else:
        # print("Please add the file you want to Encode as an Argument")
    except Exception as e:
        raise e
    # input("Press Enter to continue ...")
