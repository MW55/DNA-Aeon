#!/usr/bin/python
# -*- coding: latin-1 -*-
import math
import time
import numpy as np
from random import random

from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.OnlineBPDecoder import OnlineBPDecoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.LTBPDecoder import LTBPDecoder
from norec4dna.LTDecoder import LTDecoder
from norec4dna.distributions.IdealSolitonDistribution import IdealSolitonDistribution
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.OnlineDecoder import OnlineDecoder


def main(file):
    print("###### LT Codec number_of_chunks = 700 ######")
    print("##### Ideal Soliton Distribution #####")
    print("### LT without prioritized Packets ###")
    start = time.time()
    number_of_chunks = 700
    dist = IdealSolitonDistribution(number_of_chunks, seed=2)
    encoder = LTEncoder(file, number_of_chunks, dist)
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### LT without prioritized Packets with pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 700
    dist = IdealSolitonDistribution(number_of_chunks, seed=2)
    pseudo = LTBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(file, number_of_chunks, dist, pseudo_decoder=pseudo)
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### LT with prioritized Packets ###")
    start = time.time()
    number_of_chunks = 700
    dist = IdealSolitonDistribution(number_of_chunks, seed=2)
    encoder = LTEncoder(file, number_of_chunks, dist, prioritized_packets=[1, 3, 5])
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### LT with prioritized Packets AND pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 700
    dist = IdealSolitonDistribution(number_of_chunks, seed=2)
    pseudo = LTBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo, prioritized_packets=[1, 3, 5]
    )
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )
    print("### LT with prioritized Packets AND Gauss-pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 700
    dist = IdealSolitonDistribution(number_of_chunks, seed=2)
    pseudo = LTDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo, prioritized_packets=[1, 3, 5]
    )
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("##### Robust Soliton Distribution #####")
    print("### LT without prioritized Packets ###")
    start = time.time()
    number_of_chunks = 700
    dist = RobustSolitonDistribution(number_of_chunks, seed=2)
    encoder = LTEncoder(file, number_of_chunks, dist)
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### LT without prioritized Packets with Gauss-pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 700
    dist = RobustSolitonDistribution(number_of_chunks, seed=2)
    pseudo = LTDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(file, number_of_chunks, dist, pseudo_decoder=pseudo)
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### LT with prioritized Packets ###")
    start = time.time()
    number_of_chunks = 700
    dist = RobustSolitonDistribution(number_of_chunks, seed=2)
    encoder = LTEncoder(file, number_of_chunks, dist, prioritized_packets=[1, 3, 5])
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### LT with prioritized Packets AND Gauss-pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 700
    dist = RobustSolitonDistribution(number_of_chunks, seed=2)
    pseudo = LTDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo, prioritized_packets=[1, 3, 5]
    )
    encoder.encode_to_packets()
    end = time.time() - start
    print("Finished encoding after " + str(round(end, 4)) + " sec. " + str(
        len(encoder.get_encoded_packets())) + " Packets encoded.\n")
    print("###### Online Codec - number_of_chunks = 2500, eps=0.01, quality=3 ######")
    print("### Online without pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 2500
    epsilon = 0.01
    quality = 3
    dist = OnlineDistribution(epsilon)
    # infer number_of_chunks form Distribution:
    number_of_chunks = dist.get_size()
    pseudo = OnlineBPDecoder.pseudo_decoder()
    encoder = OnlineEncoder(file, number_of_chunks, dist, epsilon, quality)
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

    print("### Online with pseudo_decoder ###")
    start = time.time()
    number_of_chunks = 2500
    epsilon = 0.01
    quality = 3
    dist = OnlineDistribution(epsilon)
    # infer number_of_chunks form Distribution:
    number_of_chunks = dist.get_size()
    pseudo = OnlineBPDecoder.pseudo_decoder()
    encoder = OnlineEncoder(
        file, number_of_chunks, dist, epsilon, quality, pseudo_decoder=pseudo
    )
    encoder.encode_to_packets()
    end = time.time() - start
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )


# encoder -> blackbox with a (droprate) change of dropping packet -> decoder


def blackbox(encoder, decoder, droprate=0.02):
    encoder.encode_to_packets()
    i = 0
    print("[+] Created " + str(len(encoder.get_encoded_packets())) + " Packets")
    """ drop the first i elements from the set (since the set is unordered)"""
    i = math.floor(len(encoder.get_encoded_packets()) * droprate)
    packets = list(encoder.get_encoded_packets())[i:]
    for packet in packets:
        decoder.input_new_packet(packet)
    """ # OR choose fairly between all packets:    
    for packet in encoder.get_encoded_packets():
        if np.random.rand() > droprate:
            decoder.input_new_packet(packet)
        else:
            i += 1
    """
    decoder.solve()
    print("[!] Blackbox dropped " + str(i) + " random Packets.")
    return (
        decoder.is_decoded(),
        len(encoder.get_encoded_packets()),
        i,
        decoder.getSolvedCount(),
    )


def blackboxOnlineTest(file, number_of_chunks=800, droprate=0.02, seed=2, overhead=0.05):
    start = time.time()
    epsilon = 0.024912
    quality = 7
    dist = OnlineDistribution(epsilon, seed)
    number_of_chunks = dist.get_size()
    print(
        "#### Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks and a droprate of "
        + str(droprate)
        + " ####"
    )
    # pseudo = OnlineBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = OnlineEncoder(
        file, number_of_chunks, dist, epsilon, quality
    )  # , pseudo_decoder=pseudo)
    decoder = OnlineDecoder.pseudo_decoder(number_of_chunks, read_all_before_decode=True)
    encoder.set_overhead_limit(overhead)
    result, numberOfEncodedPackets, dropedCount, solvedCount = blackbox(
        encoder, decoder, droprate=droprate
    )
    end = time.time() - start
    print(
        "Blackbox-Decode "
        + ("successful" if result else "NOT successful")
        + " after "
        + str(round(end, 4))
        + " sec."
    )
    return [
        "Online_eps=" + str(epsilon) + "_quality=" + str(quality),
        result,
        numberOfEncodedPackets,
        dropedCount,
        solvedCount,
        round(end, 4),
        number_of_chunks,
    ]


def blackboxTest(file, number_of_chunks=800, droprate=0.02, seed=2, overhead=0.05):
    print(
        "#### Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks and a droprate of "
        + str(droprate)
        + " ####"
    )
    start = time.time()
    dist = RobustSolitonDistribution(S=number_of_chunks, seed=seed)
    # dist = IdealSolitonDistribution(S=number_of_chunks, seed=seed)

    # pseudo = LTBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(file, number_of_chunks, dist)  # , pseudo_decoder=pseudo)
    decoder = LTDecoder.pseudo_decoder(number_of_chunks, read_all_before_decode=True)
    encoder.set_overhead_limit(overhead)
    result, numberOfEncodedPackets, dropedCount, solvedCount = blackbox(
        encoder, decoder, droprate=droprate
    )
    end = time.time() - start
    print(
        "Blackbox-Decode "
        + ("successful" if result else "NOT successful")
        + " after "
        + str(round(end, 4))
        + " sec."
    )
    return [
        dist.get_config_string(),
        result,
        numberOfEncodedPackets,
        dropedCount,
        solvedCount,
        round(end, 4),
        number_of_chunks,
    ]


def blackboxRU10Test(
        file, number_of_chunks=800, droprate=0.02, seed=2, chunk_size=200, overhead=0.05
):
    print("Starting Blackbox Test with " + str(number_of_chunks) + " Chunks")
    start = time.time()
    dist = RaptorDistribution(number_of_chunks)
    encoder = RU10Encoder(file, number_of_chunks, dist)  # , chunk_size=chunk_size)
    decoder = RU10Decoder.pseudo_decoder(number_of_chunks, read_all_before_decode=True)
    encoder.set_overhead_limit(overhead)
    result, numberOfEncodedPackets, droped_count, solved_count = blackbox(
        encoder, decoder, droprate=droprate
    )
    end = time.time() - start
    print(
        "Blackbox-Decode "
        + ("successful" if result else "NOT successful")
        + " after "
        + str(round(end, 4))
        + " sec."
    )
    return [
        dist.get_config_string(),
        result,
        numberOfEncodedPackets,
        droped_count,
        solved_count,
        round(end, 4),
        number_of_chunks,
    ]


def get_random_int(maxInt):
    return int(random() * maxInt)


if __name__ == "__main__":
    csv = [
        "filename, codecName, number_of_chunks, numberOfEncodedPackets, droprate, dropedCount, seed, result, timeNeeded"
    ]
    file = "../.INFILES/b_lq.webm"  # "b_mq.webm"
    # main(file)
    name = "ERROR"
    for droprate in np.arange(0.04, 0.051, 0.001):
        for repeat in range(5):
            # try:
            rnd = get_random_int(math.pow(2, 31) - 1)
            number_of_chunks = 700
            name, result, numberOfEncodedPackets, dropedCount, solvedCount, timeNeeded, number_of_chunks = blackboxRU10Test(
                file, number_of_chunks=number_of_chunks, droprate=droprate, seed=rnd
            )
            line = (
                    str(file)
                    + ", "
                    + str(name)
                    + ", "
                    + str(number_of_chunks)
                    + ", "
                    + str(numberOfEncodedPackets)
                    + ", "
                    + str(droprate)
                    + ", "
                    + str(dropedCount)
                    + ", "
                    + str(solvedCount)
                    + ", "
                    + str(rnd)
                    + ", "
                    + str(result)
                    + ", "
                    + str(timeNeeded)
            )
            print(line)
            csv.append(line)
    dtimeno = "sim" + str(time.strftime("%Y-%m-%d_%H-%M", time.localtime())) + ".csv"
    with open(dtimeno, "w") as f:
        for line in csv:
            f.write(line + "\n")
else:
    print("[!] WARNING: RUN THIS SCRIPT DIRECTLY!")
