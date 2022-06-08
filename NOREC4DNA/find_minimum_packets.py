import os
import glob
import math
import bisect
import struct
import random
import argparse
import typing
from zipfile import ZipFile

import progressbar
import multiprocessing
from functools import partial
from math import ceil, floor

from norec4dna.Encoder import Encoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.Packet import ParallelPacket
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper.RepeatedTimer import RepeatedTimer
from helpful_scripts.automatedfindminimum import AutomatedFindMinimum
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode, crc32, reed_solomon_encode, dna_reed_solomon_encode
from norec4dna.helper import should_drop_packet, split_file, number_to_base_str, find_ceil_power_of_four, \
    merge_folder_content
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution

DEFAULT_ID_LEN_FORMAT = "H"
DEFAULT_NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
DEFAULT_SAVE_AS_FASTA = False
DEFAULT_SAVE_AS_ZIP = False
ONLINE_QUALITY = 7
ONLINE_EPS = 0.068

# DEFAULT_CHUNK_SIZE = 34
# NUMBER_OF_PACKETS_TO_CREATE = 655360

# global counter for progressbar
counter = None
progress_bar = None


def run(seq_seed=None, file='logo.jpg', repair_symbols=2, insert_header=False,
        error_correction=reed_solomon_encode, save_number_of_chunks_in_packet=False, l_size=1000, while_count=1000,
        chunk_size=0, number_of_chunks=300, prepend="", append="", seed_len_format=DEFAULT_ID_LEN_FORMAT,
        number_of_chunks_len_format=DEFAULT_NUMBER_OF_CHUNKS_LEN_FORMAT, method='RU10',
        mode1bmp=False, drop_above=0.4, packets_to_create=None):
    global counter
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    dna_rules = FastDNARules()
    if packets_to_create is None:
        packets_to_create = math.pow(2, 8 * struct.calcsize(seed_len_format))
    rules = dna_rules
    if repair_symbols != 0:
        dist, error_correction = get_err_dist(method, number_of_chunks, repair_symbols)
    else:
        dist = RaptorDistribution(number_of_chunks)
    if method == 'RU10':
        x = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=insert_header, rules=rules,
                        error_correction=error_correction, id_len_format=seed_len_format,
                        number_of_chunks_len_format=number_of_chunks_len_format,
                        save_number_of_chunks_in_packet=save_number_of_chunks_in_packet, mode_1_bmp=mode1bmp,
                        prepend=prepend, append=append)
        x.prepare()
    elif method == 'LT':
        x = LTEncoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=insert_header, rules=rules,
                      error_correction=error_correction, number_of_chunks_len_format=number_of_chunks_len_format,
                      id_len_format=seed_len_format, save_number_of_chunks_in_packet=save_number_of_chunks_in_packet)
        x.prepareEncoder()
    elif method == 'Online':
        number_of_chunks = dist.get_size()
        x = OnlineEncoder(file, number_of_chunks, dist, ONLINE_EPS, ONLINE_QUALITY, error_correction=error_correction,
                          quality_len_format="B", insert_header=False, check_block_number_len_format=seed_len_format,
                          number_of_chunks_len_format=number_of_chunks_len_format, rules=rules,
                          save_number_of_chunks_in_packet=False)
        x.prepare()
    else:
        raise NotImplementedError("Choose: RU10, LT or Online")
    i = 0
    tmp_list = []
    while i < while_count:
        if seq_seed is not None:
            if seq_seed + i >= packets_to_create:
                break
            packet = x.create_new_packet(seed=seq_seed + i)
        else:
            packet = x.create_new_packet()
        if i == 0:
            print(f"%i , %s" % (len(packet.get_dna_struct(True)), packet.get_dna_struct(True)))
        _ = should_drop_packet(rules, packet)
        if packet.error_prob <= drop_above and (len(tmp_list) < l_size or packet.error_prob < tmp_list[-1].error_prob):
            if packet not in tmp_list:
                bisect.insort_left(tmp_list, packet)
            else:
                elem = next((x for x in tmp_list if x == packet), None)
                if packet < elem:
                    tmp_list.remove(elem)
                    del elem
                    bisect.insort_left(tmp_list, packet)
            if len(tmp_list) > l_size:
                for ele1m in tmp_list[l_size + 1:]:
                    del ele1m
                tmp_list = tmp_list[:l_size]

        else:
            del packet
        i += 1
        # += operation is not atomic, so we need to get a lock:
        with counter.get_lock():
            counter.value += 1
    # save_packets_fasta(tmp_list, out_file=method + "_out_partial", file_ending="." + method + "_DNA",
    #                   clear_output=False)
    conf = {'error_correction': error_correction, 'repair_symbols': repair_symbols,
            'number_of_splits': _number_of_splits,
            'find_minimum_mode': True, 'seq_seed': seq_seed}
    # x.save_config_file(conf, section_name=method + "_" + file)
    if x.progress_bar is not None:
        x.progress_bar.finish()
    return [ParallelPacket.from_packet(p) for p in tmp_list]


def save_packets_zip(encodedPackets, out_file: typing.Optional[str] = None, file_ending=".zip", seed_is_filename=True):
    if not out_file.endswith(file_ending):
        out_file = out_file + "." + file_ending
    i = 0
    abs_dir = os.path.split(os.path.abspath("../" + out_file))[0]
    if not os.path.exists(abs_dir):
        os.makedirs(abs_dir)

    with ZipFile(out_file, 'w') as f:
        for packet in encodedPackets:
            if seed_is_filename:
                i = packet.id
                f.writestr(f"{i}{file_ending}", packet.get_struct(True))
            i += 1
    print(f"Saved result at: %s" % out_file)


def save_packets_fasta(packets, out_file, file_ending, clear_output=True, seed_is_filename=True):
    if not out_file.endswith("/"):
        files = glob.glob(out_file + "/*")
    else:
        files = glob.glob(out_file + "*")
    if clear_output:
        for f in files:
            try:
                os.remove(f)
            except Exception as ex:
                print("Error: ", ex)
    i = 0
    if not os.path.exists(out_file):
        os.makedirs(out_file)
    out_str = out_file + "/" + str(random.randint(1, 80000000000)) + "_" + str(
        random.randint(1, 800000000000)) + ".fasta"
    with open(out_str, "w") as f:
        for packet in packets:  # sorted(packets, key=lambda elem: (elem.error_prob, elem.__hash__())):
            if seed_is_filename:
                i = packet.id
            e_prob = (str(ceil(packet.error_prob * 100)) + "_") if packet.error_prob is not None else ""
            f.write(">" + e_prob + str(i) + file_ending + "\n" + packet.get_dna_struct(True) + "\n")
            i += 1
    print(f"Saved result at: %s" % out_str)


def save_packets(packets, out_file, file_ending, clear_output=True, seed_is_filename=True):
    if not out_file.endswith("/"):
        files = glob.glob(out_file + "/*")
    else:
        files = glob.glob(out_file + "*")
    if clear_output:
        for f in files:
            try:
                os.remove(f)
            except Exception as ex:
                print("Error: ", ex)
    i = 0
    if not os.path.exists(out_file):
        os.makedirs(out_file)
    for packet in packets:  # sorted(packets, key=lambda elem: (elem.error_prob, elem.__hash__())):
        if seed_is_filename:
            i = packet.id
        e_prob = (str(ceil(packet.error_prob * 100)) + "_") if packet.error_prob is not None else ""
        with open(out_file + "/" + e_prob + str(i) + file_ending, "w") as f:
            f.write(packet.get_dna_struct(True))
        i += 1
    print("Saved files in folder: %s" % out_file)


def reduce_lists(base, input_list=None, l_size=100):
    if input_list is None and len(base) == 2:
        input_list = base[1]
        base = base[0]
    base = base[:l_size]
    for packet in input_list:
        if packet.error_prob < base[-1].error_prob:
            if packet not in base:
                bisect.insort_left(base, packet)
            else:
                elem = next((x for x in base if x == packet), None)
                if packet < elem:
                    base.remove(elem)
                    bisect.insort_left(base, packet)
            if len(base) > l_size:
                base = base[:l_size]
        else:
            # we can skip the rest of the current result-row since all remaining packets have a higher error_prob
            # than the worst packet of the current answer
            return base
    return base


def init_mp(args):
    global counter
    counter = args


def create_progress_bar(max_value):
    widgets = [
        progressbar.Percentage(),
        progressbar.Bar(),
        ' Correct: ',
        progressbar.Counter(),
        ', ',
        progressbar.AdaptiveETA(), ' ',
        progressbar.Timer()
    ]
    return progressbar.ProgressBar(max_value=max_value, widgets=widgets, max_error=False,
                                   redirect_stdout=True).start()


def update_progressbar():
    progress_bar.update(counter.value)


# multiprocess and merge all created lists at the end
def main(filename="logo.jpg", repair_symbols=2, while_count=1000, out_size=1000, chunk_size=0, number_of_chunks=300,
         sequential=False, spare1core=False, prepend="", append="", insert_header=False,
         seed_len_format=DEFAULT_ID_LEN_FORMAT,
         method='RU10', mode1bmp=False, drop_above=0.4, save_as_fasta=DEFAULT_SAVE_AS_FASTA,
         packets_to_create=None, save_as_zip=DEFAULT_SAVE_AS_ZIP):
    global progress_bar, counter
    if packets_to_create is None:
        packets_to_create = math.pow(2, 8 * struct.calcsize(seed_len_format))
    cores = multiprocessing.cpu_count()
    if spare1core:
        cores = cores - 1
    counter = multiprocessing.Value('i', 0)
    p = multiprocessing.Pool(cores, initializer=init_mp, initargs=(counter,))
    param = [None] * cores
    if sequential:
        stepsize = packets_to_create / cores
        param = [int(floor(i * stepsize)) for i in range(cores)]
        while_count = int(ceil(stepsize)) + 1

    progress_bar = create_progress_bar(packets_to_create)
    progress_bar.update(0)
    rt = RepeatedTimer(2, update_progressbar, )
    a = p.map(
        partial(run, file=filename, repair_symbols=repair_symbols, l_size=out_size, while_count=while_count,
                chunk_size=chunk_size, number_of_chunks=number_of_chunks, prepend=prepend, append=append,
                insert_header=insert_header, seed_len_format=_seed_size_str, method=method, mode1bmp=mode1bmp,
                drop_above=drop_above), param)
    rt.stop()
    progress_bar.finish()
    print("Merging results...")
    progress_bar.start(len(a))
    progress_bar.update(0)
    while len(a) > 1:
        if len(a) % 2 != 0:
            a.append([])
        tmp = []
        for i in range(0, len(a) - 1, 2):
            tmp.append([a[i], a[i + 1]])
        a = p.map(partial(reduce_lists, input_list=None, l_size=out_size), tmp)
        progress_bar.update(progress_bar.max_value - len(a))
    a = a[0]
    rt.stop()

    out_file = "parallel_" + method
    file_ending = "_joined." + method + "_DNA"
    if save_as_fasta:
        save_packets_fasta(a, out_file=out_file, file_ending=file_ending)
    elif save_as_zip:
        save_packets_zip(a, out_file=out_file, file_ending=file_ending)
    else:
        save_packets(a, out_file=out_file, file_ending=file_ending)
    progress_bar.finish()
    return a, out_file, file_ending


def get_err_dist(_method, _number_of_chunks, _repair_symbols):
    if _method == 'RU10':
        dist = RaptorDistribution(_number_of_chunks)
    elif _method == 'LT':
        dist = ErlichZielinskiRobustSolitonDistribution(_number_of_chunks, seed=2)
    elif _method == 'Online':
        dist = OnlineDistribution(ONLINE_EPS)
    else:
        raise NotImplementedError("Choose: RU10, LT or Online")
    return dist, lambda x: reed_solomon_encode(x, _repair_symbols)


def init_optimization(packets, error_correction, dist):
    return AutomatedFindMinimum(packets, error_correction=error_correction, dist=dist)


def instantiate_error_correction(err_corr_name, no_repair_symbols):
    if err_corr_name == "nocode":
        return nocode
    elif err_corr_name == "crc":
        return crc32
    elif err_corr_name == "reedsolomon":
        if no_repair_symbols != 2:
            return lambda x: reed_solomon_encode(x, no_repair_symbols)
        else:
            return reed_solomon_encode
    elif err_corr_name == "dna_reedsolomon":
        if no_repair_symbols != 2:
            return lambda x: dna_reed_solomon_encode(x, no_repair_symbols)
        else:
            return dna_reed_solomon_encode
    else:
        raise RuntimeError(
            "Selected Error Correction not supported, choose: 'nocode', 'crc', 'reedsolomon' or 'dna_reedsolomon")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--list_size", metavar="list_size", required=False, type=int, default=1000,
                        help="size of operational list per thread [inferred by #cores if sequential is set to true]")
    parser.add_argument("--out_size", metavar="out_size", required=False, type=int, default=1000,
                        help="number of packets to save")
    parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=0,
                        help="chunk size (default=0 [infer from number of chunks])")
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=False, type=int, default=300,
                        help="number of chunks (default=300) [ignored if chunk_size is set to value != 0]")
    parser.add_argument("--sequential", required=False, default=False, action="store_true")
    parser.add_argument("--spare1core", required=False, default=False, action="store_true")
    parser.add_argument("--split_input", metavar="split_input", required=False, type=int, default=1,
                        help="number of sub-codes to split input file into")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--insert_header", required=False, default=False, action="store_true")
    parser.add_argument("--method", required=False, type=str, default="RU10",
                        help="RU10, Online or LT are supported")
    parser.add_argument("--optimization", required=False, default=False, action="store_true")
    parser.add_argument("--overhead", required=False, type=float, default=0.2,
                        help="overhead to use for the optimization")
    parser.add_argument("--overhead_factor", required=False, type=float, default=0.0,
                        help="factor to exceed the overhead with if necessary to optimize")
    parser.add_argument("--error_prob_factor", required=False, type=float, default=0.1,
                        help="factor for the max error_prob for new packets based on the avg of the existing ones")
    parser.add_argument("--plot", required=False, default=False, action="store_true",
                        help="plot the results of the optimization")
    parser.add_argument("--mode1bmp", required=False, default=False, action="store_true",
                        help="input is an image, save as b/w bmp without header (applies to RU10 only")
    parser.add_argument("--drop_above", required=False, type=float, default=1.0,
                        help="Instantly drop packet if error-prob above set threshold")
    parser.add_argument("--seed_size_str", required=False, type=str, default=DEFAULT_ID_LEN_FORMAT,
                        help="struct-string for seed - possible values: I,H,B (see struct) This impacts the amount of generated packets on sequential mode")
    parser.add_argument("--store_as_fasta", required=False, default=DEFAULT_SAVE_AS_FASTA,
                        action="store_true", help="if set, store result as fasta file")
    parser.add_argument("--store_as_zip", required=False, default=DEFAULT_SAVE_AS_ZIP,
                        action="store_true", help="if set, store result as zip file")
    """
    parser.add_argument("--savenumberofchunks", metavar="savenumberofchunks", required=False, type=bool,
                        default=False)
    """
    args = parser.parse_args()
    _input_file = args.filename
    _repair_symbols = args.repair_symbols
    _list_size = args.list_size
    _out_size = args.out_size
    _chunk_size = args.chunk_size
    _number_of_hunks = args.number_of_chunks
    _sequential = args.sequential
    _spare1core = args.spare1core
    _number_of_splits = args.split_input
    _insert_header = args.insert_header
    _error_correction = args.error_correction
    _drop_above = args.drop_above
    _seed_size_str = args.seed_size_str
    _store_as_fasta = args.store_as_fasta
    _store_as_zip = args.store_as_zip
    _method = args.method
    _mode1bmp = args.mode1bmp

    _optimization = args.optimization
    _overhead = args.overhead
    _overhead_factor = args.overhead_factor
    _error_prob_factor = args.error_prob_factor
    _plot = args.plot

    _error_correction = instantiate_error_correction(_error_correction, _repair_symbols)

    if _chunk_size != 0:
        _number_of_hunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(_input_file, _chunk_size)
    if _list_size < _out_size:
        print("Warning: A list size smaller than the out_size may result in non-optimal results.")
    if _optimization:
        print("Warning: Please make sure to generate enough packets for the automated optimization.")

    if not _optimization:
        if _number_of_splits > 1:
            input_files = split_file(_input_file, _number_of_splits)
            power_of_four = find_ceil_power_of_four(len(input_files))
            print("Spltting input into {} sub files. We need to prepend/append {} base(s)".format(len(input_files),
                                                                                                  power_of_four))
            prepend_matching = {input_files[i]: number_to_base_str(i, power_of_four) for i in range(len(input_files))}
        else:
            input_files = [_input_file]
            prepend_matching = {_input_file: ""}
        for _input_file in input_files:
            try:
                print("File to encode: " + str(_input_file))
                #def callab():
                main(_input_file, _repair_symbols, _list_size, _out_size, _chunk_size, _number_of_hunks,
                         _sequential, _spare1core, method=_method, append=prepend_matching[_input_file],
                         insert_header=_insert_header, seed_len_format=_seed_size_str, drop_above=_drop_above,
                         save_as_fasta=_store_as_fasta, save_as_zip=_store_as_zip)


                #print(timeit.repeat(callab, number=1111, repeat=5))
            except Exception as ex:
                print(ex)
                raise ex
                #list_fds()
        # if len(input_files) > 1:
        #    merge_folder_content("split_" + os.path.basename(_input_file), _input_file + "_combined_split_output",
        #                         append_folder_name=True,
        #                         clear_dest_folder=True)
    else:
        par_packets, _out_file, _file_ending = main(_input_file, _repair_symbols, _list_size, _out_size, _chunk_size,
                                                    _number_of_hunks, _sequential, _spare1core, method=_method,
                                                    seed_len_format=_seed_size_str, drop_above=_drop_above)
        distribution, _error_correction = get_err_dist(_method, _number_of_hunks, _repair_symbols)
        optimizer = init_optimization(par_packets, _error_correction, distribution)
        optimized_packets = optimizer.automated_optimization(_overhead, _overhead_factor, _error_prob_factor,
                                                             plot=_plot)
        save_packets(optimized_packets, out_file=_out_file + "optimized/", file_ending=_file_ending)
