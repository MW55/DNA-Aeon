import os
import time
import math
import struct
import numpy as np
from random import random

from norec4dna.helper import xor_mask
from norec4dna.rules.DNARules import DNARules
from norec4dna.helper.bin2Quaternary import string2QUATS
from norec4dna.helper.quaternary2Bin import quats_to_bytes

lines = [
    "Algorithm,A_Permutation,T_Permutation,C_Permutation,G_Permutation,dinucleotid_Runs,Homopolymers,GC_Content,Trinucleotid_Runs,Random_Permutation,Overall_Dropchance,Random_Number,Did_Drop"
]


def get_file_size(file):
    return os.stat(file).st_size


def get_number_of_chunks_for_file_with_chunk_size(file, chunk_size, insert_header=False):
    file_size = get_file_size(file)
    return math.ceil(1.0 * file_size / chunk_size) + 1 if insert_header else 0


def get_random_int(max_int):
    return int(random() * max_int)


def create_packet(leng=100):
    tmp = np.random.bytes(
        leng + 16
    )  # +16 to adapt for average header in raptor/online/lt...
    return tmp


def should_drop_packet(packet, add_line=True, scale=1.0):
    rand = random()
    dna_data = "".join(string2QUATS(packet))
    drop_chance, data, _ = DNARules.apply_all_rules_with_data(dna_data)
    drop_chance = scale * drop_chance
    if add_line:
        line = ("random_bytes," + ",".join([str(round(x, 4)) for x in data]) + "," + str(drop_chance) + "," + str(
            rand) + "," + str(drop_chance > rand))
        lines.append(line)
    # drop packet if rand bigger than the drop_chance for this Packet.
    return drop_chance > rand


def blackbox(file, number_of_chunks, seed, overhead, leng=100):
    start = time.time()
    name = "random_bytes"
    no_packets = int(math.ceil(number_of_chunks * (1.0 + overhead)))
    packets = [create_packet(leng) for _ in range(no_packets)]
    result = True
    correct = 0
    for elem in packets:
        elem = struct.pack(
            "<I" + str(len(elem)) + "sI", xor_mask(leng), elem, xor_mask(leng)
        )  # simulates a crc
        if should_drop_packet(elem, True):
            result = False
        else:
            correct += 1
    end = time.time() - start

    return name, result, number_of_chunks, number_of_chunks - correct, round(end, 4)


def main(file="logo.jpg", repeats=50):
    csv = [
        "filename, overhead, codecName, number_of_chunks, invalid_drop, seed, result, time_needed"
    ]
    # main(file)
    name = "ERROR"
    mode = "random_bytes"
    file_size = get_file_size(file)
    for overhead in [0.0]:  # np.arange(0.05, 0.1,
        #     0.005):  # [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]:

        chunk_size = 100
        insert_header = True
        if chunk_size != 0:
            number_of_chunks = get_number_of_chunks_for_file_with_chunk_size(
                file, chunk_size=chunk_size, insert_header=insert_header
            )
        else:
            number_of_chunks = 800
        number_repair_symbols = int(math.floor((file_size // number_of_chunks) * overhead))
        for _ in range(repeats):
            rnd = get_random_int(math.pow(2, 31) - 1)

            try:
                name, result, number_of_chunks, invalid_drop, time_needed = blackbox(
                    file,
                    number_of_chunks=number_of_chunks,
                    seed=rnd,
                    overhead=overhead,
                    leng=chunk_size,
                )
                line = (str(file) + "," + str(overhead) + "," + str(name) + "," + str(number_of_chunks) + "," + str(
                    invalid_drop) + "," + str(rnd) + "," + str(result) + "," + str(time_needed))
            except Exception:
                line = (str(file) + "," + str(overhead) + "," + str(name) + "," + str(
                    number_of_chunks) + "," + "ERROR" + "," + str(rnd) + "," + "ERROR" + "," + "ERROR")
            print(line)
            csv.append(line)

        dtimeno = (mode + "_" + str(number_repair_symbols) + "_" + str(overhead) + "_sim" + str(
            time.strftime("%Y-%m-%d_%H-%M", time.localtime())) + ".csv")

        with open("DNA_" + dtimeno, "w") as f:
            for line in lines:
                f.write(line + "\n")
        lines.clear()
        lines.append(
            "Algorithm,A_Permutation,T_Permutation,C_Permutation,G_Permutation,dinucleotid_Runs,Homopolymers,"
            "GC_Content,Trinucleotid_Runs,Random_Permutation,Overall_Dropchance,Random_Number,Did_Drop"
        )

        with open(dtimeno, "w") as f:
            for line in csv:
                f.write(line + "\n")
        csv = [
            "filename, overhead, codecName, number_of_chunks, invalid_drop, seed, result, time_needed"
        ]


if __name__ == "__main__":
    main()
