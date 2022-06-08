#!/usr/bin/python
import argparse
from norec4dna.LTDecoder import LTDecoder
from norec4dna.ErrorCorrection import nocode, get_error_correction_decode
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution

STATIC_NUM_CHUNKS = 31
NULL_IS_TERMINATOR = False
IMPLICIT_MODE = True
HEADER_CHUNK = True
PRINT_TO_OUTPUT = False
NUMBER_OF_CHUNKS_LEN_FORMAT = "H"
ID_LEN_FORMAT = "H"
CRC_LEN_FORMAT = "L"
PACKET_LEN_FORMAT = "I"


# NUMBER_OF_CHUNKS_IN_PACKET := STATIC_NUM_CHUNKS is None

# --error_correction reedsolomon --repairsymbols 2 LT_Dorn.tar.gz (num_chunks == 69)
# --error_correction reedsolomon --repairsymbols 2 LT_Dorn (num_chunks == 153)
class demo_decode:
    @staticmethod
    def decode(file, error_correction=nocode, null_is_terminator=False, mode_1_bmp=False,
               number_of_chunks=STATIC_NUM_CHUNKS, use_header_chunk=False, id_len_format=ID_LEN_FORMAT,
               number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT, packet_len_format=PACKET_LEN_FORMAT,
               crc_len_format=CRC_LEN_FORMAT):
        dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)
        """try:
            decoder = LTBPDecoder(file, error_correction=error_correction, use_headerchunk=HEADER_CHUNK,
                                      static_number_of_chunks=STATIC_NUM_CHUNKS, implicit_mode=IMPLICIT_MODE, dist=dist)
            print("[1/2] Approximation Decode")
            _internal(decoder)
            if not decoder.is_decoded():
                print("[2/2] Approximation Decode")
                _internal(decoder)
            decoder.saveDecodedFile(null_is_terminator=NULL_IS_TERMINATOR)
        except Exception as e:"""
        print("[X/2] Falling back to Gauss-Mode")
        print("Falling back to Gauss-Mode")
        decoder = LTDecoder(file, error_correction=error_correction, use_headerchunk=use_header_chunk,
                            static_number_of_chunks=number_of_chunks, implicit_mode=IMPLICIT_MODE, dist=dist)
        decoder.read_all_before_decode = True
        decoder.decode(number_of_chunks_len_format=number_of_chunks_len_format, seed_len_format=id_len_format,
                       degree_len_format="H")
        decoder.solve()
        decoder.saveDecodedFile(null_is_terminator=null_is_terminator, print_to_output=PRINT_TO_OUTPUT)
        print(decoder.tmp_mapping)


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("filename", metavar="file", type=str, help="the file / folder to Decode")
        parser.add_argument("--error_correction", metavar="error_correction", type=str, required=False,
                            default="nocode",
                            help="Error Correction Method to use; possible values: \
                                nocode, crc, reedsolomon (default=nocode)")
        parser.add_argument("--repair_symbols", metavar="repair_symbols", type=int, required=False, default=2,
                            help="number of repair symbols for ReedSolomon (default=2)")
        args = parser.parse_args()
        filename = args.filename
        e_correction_str = args.error_correction
        norepair_symbols = args.repair_symbols
        error_correction = get_error_correction_decode(e_correction_str, norepair_symbols)
        print("Zu dekodierende Datei / Ordner: " + str(filename))
        demo = demo_decode()
        demo.decode(filename, error_correction=error_correction)
    except Exception as e:
        print(e)
    # input("Press Enter to continue ...")
