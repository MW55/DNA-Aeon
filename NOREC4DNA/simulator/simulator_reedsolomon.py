import time
import math
import struct
import numpy as np
from random import random

from norec4dna.rules.DNARules import DNARules
from norec4dna.helper.bin2Quaternary import string2QUATS
from norec4dna.helper.quaternary2Bin import quats_to_bytes
from norec4dna.ReedSolomonSuite import get_file_size, ReedSolomonEncoder, ReedSolomonDecoder, xor_mask

lines = [
    "Algorithm,CRC,A_Permutation,T_Permutation,C_Permutation,G_Permutation,dinucleotid_Runs,Homopolymers,GC_Content,Trinucleotid_Runs,Random_Permutation,Overall_Dropchance,Random_Number,Did_Drop"
]


def get_number_of_chunks_for_file_with_chunk_size(file, chunk_size, insert_header=False):
    file_size = get_file_size(file)
    return math.ceil(1.0 * file_size / chunk_size) + 1 if insert_header else 0


def get_random_int(max_int):
    return int(random() * max_int)


def permute(data, add_line=True):
    dna_data = "".join(string2QUATS(data))
    drop_chance, data, _ = DNARules.apply_all_rules_with_data(dna_data)
    drop_chance = 0.01 * drop_chance
    if add_line:
        line = ("ReedSolomon" + "," + str() + "," + ",".join([str(round(x, 4)) for x in data]) + "," + str(drop_chance))
        lines.append(line)
    lim_drop_chance = min(drop_chance, 0.9999999)
    number_of_permutations = int(math.ceil(lim_drop_chance * len(dna_data)))
    dna_data = list(dna_data)
    for i in range(number_of_permutations):
        dna_data[i] = flip_base(dna_data[i])  # flip the dna-base to an incorrrect one
    data = b""
    for x in list(zip(*[iter(dna_data)] * 4)):
        data += quats_to_bytes("".join(x))
    return data


def flip_base(base):
    if base == "A":
        return "C"
    elif base == "C":
        return "A"
    elif base == "G":
        return "T"
    else:
        return "G"


def blackbox(file, number_of_chunks, seed, overhead, r_symbols):
    start = time.time()
    name = "ReedSolomon"
    rs = ReedSolomonEncoder(file, number_of_chunks, 0, overhead)
    a = rs.create_chunks()

    dec = ReedSolomonDecoder()
    resmap = {}
    result = True
    correct = 0
    for elem in a:
        # permute input according to dna rules:
        elem = permute(elem)
        elem = struct.pack(
            "<I" + str(len(elem)) + "s", xor_mask(rs.number_repair_symbols), elem
        )
        try:
            i, res = dec.decode(elem)
            resmap[i] = res
            correct += 1
        except:
            result = False
    end = time.time() - start
    return name, result, number_of_chunks, number_of_chunks - correct, round(end, 4)


def main(file="logo.jpg", repeats=5):
    csv = [
        "filename, overhead, codecName, number_of_chunks, invalid_drop, seed, result, time_needed"
    ]
    name = "ERROR"
    mode = "ReedSolomon"
    file_size = get_file_size(file)
    for overhead in np.arange(0.05, 0.1, 0.005):
        chunk_size = 100
        insert_header = True
        if chunk_size != 0:
            number_of_chunks = get_number_of_chunks_for_file_with_chunk_size(
                file, chunk_size=chunk_size, insert_header=insert_header
            )
        else:
            number_of_chunks = 800
        number_repair_symbols = int(math.floor((file_size // number_of_chunks) * overhead))
        print("No_repair_symbols=" + str(number_repair_symbols))
        for _ in range(repeats):
            rnd = get_random_int(math.pow(2, 31) - 1)

            try:
                name, result, number_of_chunks, invalid_drop, time_needed = blackbox(
                    file,
                    number_of_chunks=number_of_chunks,
                    seed=rnd,
                    overhead=overhead,
                    r_symbols=number_repair_symbols,
                )
                line = (
                        str(file)
                        + ","
                        + str(overhead)
                        + ","
                        + str(name)
                        + ","
                        + str(number_of_chunks)
                        + ","
                        + str(invalid_drop)
                        + ","
                        + str(rnd)
                        + ","
                        + str(result)
                        + ","
                        + str(time_needed)
                )
            except Exception:
                line = (
                        str(file)
                        + ","
                        + str(overhead)
                        + ","
                        + str(name)
                        + ","
                        + str(number_of_chunks)
                        + ","
                        + "ERROR"
                        + ","
                        + str(rnd)
                        + ","
                        + "ERROR"
                        + ","
                        + "ERROR"
                )
            print(line)
            csv.append(line)

        dtimeno = (
                mode
                + "_"
                + str(number_repair_symbols)
                + "_"
                + str(overhead)
                + "_sim"
                + str(time.strftime("%Y-%m-%d_%H-%M", time.localtime()))
                + ".csv"
        )
        with open("DNA_" + dtimeno, "w") as f:
            for line in lines:
                f.write(line + "\n")
        lines.clear()
        lines.append(
            "Algorithm,CRC,A_Permutation,T_Permutation,C_Permutation,G_Permutation,dinucleotid_Runs,Homopolymers,"
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
