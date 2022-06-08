#!/usr/bin/python
import os
import glob
import argparse
import progressbar
import multiprocessing

from norec4dna.Encoder import Encoder
from norec4dna.Packet import ParallelPacket
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode

CHUNK_SIZE = 100


def create_progress_bar(max_value):
    widgets = [
        progressbar.Percentage(),
        progressbar.Bar(),
        ' Correct: ',
        progressbar.Counter(),
        ', ',
        progressbar.Variable('Corrupt'),
        ', ',
        progressbar.AdaptiveETA(), ' ',
        progressbar.Timer()
    ]
    return progressbar.ProgressBar(max_value=max_value, widgets=widgets, max_error=False,
                                   redirect_stdout=True).start()


def encode(p_output, file, as_dna=True, error_correction=nocode, insert_header=False,
           save_number_of_chunks_in_packet=False, overhead=6.0, clear_output=False):
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, CHUNK_SIZE)
    dist = RaptorDistribution(number_of_chunks)
    dna_rules = FastDNARules()
    if as_dna:
        rules = dna_rules
    else:
        rules = None
    x = RU10Encoder(file, number_of_chunks, dist, chunk_size=CHUNK_SIZE, insert_header=insert_header, rules=rules,
                    error_correction=error_correction, id_len_format="H", number_of_chunks_len_format="B",
                    save_number_of_chunks_in_packet=save_number_of_chunks_in_packet)
    x.set_overhead_limit(overhead)
    x.encode_to_packets()
    p_output.send([ParallelPacket.from_packet(packet) for packet in x.encodedPackets])
    p_output.send("DONE")
    p_output.close()
    return 0


if __name__ == "__main__":
    # try:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--asdna",
        help="convert packets to dna and use dna rules",
        action="store_true",
        required=False,
    )
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", metavar="insert_header", required=False, type=bool, default=False)
    parser.add_argument("--save_number_of_chunks", metavar="save_number_of_chunks", required=False, type=bool,
                        default=False)
    args = parser.parse_args()
    _overhead = 6.0
    _file = args.filename
    _as_dna = args.asdna
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _save_number_of_chunks = args.save_number_of_chunks
    _error_correction = get_error_correction_encode(args.error_correction, _repair_symbols)
    print("File to encode: " + str(_file))

    from multiprocessing import Process

    cores = multiprocessing.cpu_count()

    print("[i] Clearing output-folder...")
    fulldir, filename = os.path.split(os.path.realpath(_file))
    filename = "RU10_" + filename
    out_file = os.path.join(fulldir, filename)
    if not out_file.endswith("/"):
        files = glob.glob(out_file + "/*")
    else:
        files = glob.glob(out_file + "*")
    for f in files:
        os.remove(f)

    print("[i] Spawning " + str(cores) + " processes:")
    processes = []
    for core in range(cores):
        _p_output, _p_input = multiprocessing.Pipe()
        p = Process(target=encode, args=(
            _p_output, _file, _as_dna, _error_correction, False, False,
            -(1.0 - (1.0 / cores)) + 1.0 * _overhead / cores,
            False))
        p.start()
        print("[" + str(core + 1) + "] started")
        processes.append((p, _p_input))

    res = set()
    for process, _p_input in processes:
        done = False
        while not done:
            inp = _p_input.recv()
            if inp != "DONE":
                res = res.union(inp)
            else:
                _p_input.close()
                done = True
        process.join()

    _number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(_file, CHUNK_SIZE)
    _dist = RaptorDistribution(_number_of_chunks)
    tmp = RU10Encoder(_file, _as_dna, distribution=_dist, chunk_size=CHUNK_SIZE)
    tmp.encodedPackets = res
    tmp.save_packets(True, save_as_dna=_as_dna, clear_output=True, seed_is_filename=True)

    # input("Press Enter to continue ...")
