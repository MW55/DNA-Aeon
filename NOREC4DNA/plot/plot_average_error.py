import struct
import numpy as np
import multiprocessing
from functools import partial

from norec4dna import RU10Encoder, reed_solomon_encode
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper import should_drop_packet

SEED_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
save_number_of_chunks_in_packet = True
NO_REPAIR_SYMBOLS = 4
INSERT_HEADER = True


def get_error_sum(file, number_of_chunks, chunk_size, seq_seed=None, while_count=1000):
    max_seed = np.power(2, 8 * struct.calcsize(SEED_LEN_FORMAT))
    dist = RaptorDistribution(number_of_chunks)
    dna_rules = FastDNARules()
    error_correction = lambda x: reed_solomon_encode(x, NO_REPAIR_SYMBOLS)
    encoder = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=INSERT_HEADER,
                          rules=dna_rules, error_correction=error_correction, id_len_format=SEED_LEN_FORMAT,
                          number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
                          save_number_of_chunks_in_packet=save_number_of_chunks_in_packet, prepend="", append="")
    encoder.prepare()
    i = 0
    res = []
    while i < while_count:
        if seq_seed is not None:
            if seq_seed + i >= max_seed:
                break
            packet = encoder.create_new_packet(seed=seq_seed + i)
        else:
            packet = encoder.create_new_packet()
        should_drop_packet(dna_rules, packet)
        res.append(packet.error_prob)
    return res


def main(file, number_of_chunks, chunk_size, spare1core=True, sequential=True, while_count=1000):
    cores = multiprocessing.cpu_count()
    if spare1core:
        cores = cores - 1
    p = multiprocessing.Pool(cores)
    param = [None] * cores
    if sequential:
        stepsize = 65536 / cores
        param = [int(np.floor(i * stepsize)) for i in range(cores)]
        while_count = int(np.ceil(stepsize)) + 1
    a = p.map(partial(get_error_sum, file, number_of_chunks, chunk_size, while_count=while_count), param)
    print(a)


if __name__ == "__main__":
    main("logo.jpg", 400, 0, True, True)
