import math
import multiprocessing
import struct
import sys
import typing
from functools import partial
from numpy import mean, ndarray
from scipy.optimize import minimize

from norec4dna import RU10Encoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper import should_drop_packet
from norec4dna.rules.FastDNARules import FastDNARules
from optimization_helper import list_to_diff_list, scale_to, diff_list_to_list

__FILE = "logo.jpg"
__NUM_CHUNKS = 140
__NUM_CHUNKS_LEN_FORMAT = "H"
__runs = 0


def read_sig(signal_number, frame):
    print("File saved, process continuing.")


def kill_sig(signal_number, frame):
    print("File saved, process killed.")
    sys.exit()


def create_packets_e_prob(start_num: int, normed_dist: ndarray, number_of_packets: int, rules=None):
    dist_obj = RaptorDistribution(__NUM_CHUNKS)
    dist_obj.f = normed_dist
    dist_obj.d = [x for x in range(0, 41)]
    encoder = RU10Encoder(file=__FILE, number_of_chunks=__NUM_CHUNKS, distribution=dist_obj, insert_header=False)
    encoder.prepare()
    if rules is None:
        rules = FastDNARules()
    packets_e_prob = []
    for i in range(start_num, start_num + number_of_packets):
        packet = encoder.create_new_packet(seed=i)
        should_drop_packet(rules, packet)
        packets_e_prob.append(packet.error_prob)
        del packet
    del encoder
    return packets_e_prob


def norm_list(lst):
    y = sum(lst)
    return [x / y for x in lst]


def calculate_average_error_for_distribution(distribution_list: typing.List[int]):
    global __runs
    __runs += 1
    print(__runs)
    cores = multiprocessing.cpu_count() - 1
    p = multiprocessing.Pool(cores)
    stepsize = math.floor(math.pow(2, struct.calcsize(__NUM_CHUNKS_LEN_FORMAT) * 8)) / cores
    param = [int(math.floor(i * stepsize)) for i in range(cores)]
    while_count = int(math.ceil(stepsize)) + 1
    distribution_list = scale_to(norm_list(diff_list_to_list(distribution_list)), 1048576)
    print(distribution_list)
    e_prob = p.map(
        partial(create_packets_e_prob, normed_dist=norm_list(distribution_list), number_of_packets=while_count),
        param)
    p.close()
    p.join()
    e_prob = sorted([item for sublist in e_prob for item in sublist])[:__NUM_CHUNKS * 4]
    return mean(e_prob)  # sort items and then take the first X and mean them.


def opt_callback(xk, res):
    print(xk)
    print(res)
    return False


def minimize_nm(dist, file, hyper_parameters=None):
    dist = list_to_diff_list(dist)
    res = minimize(calculate_average_error_for_distribution, dist, method='nelder-mead', callback=opt_callback)
    print(res)


if __name__ == "__main__":
    minimize_nm([0, 10241, 491582, 712794, 831695, 831695, 831695, 831695, 831695, 831695, 948446, 1032189, 1032189,
                 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189,
                 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189,
                 1032189, 1032189, 1032189, 1032189, 1032189, 1048576], "Dorn")
# implements Nelder-mead
