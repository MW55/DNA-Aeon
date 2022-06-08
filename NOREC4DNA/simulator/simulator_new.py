#!/usr/bin/python
# -*- coding: latin-1 -*-
import os
import time
import math
import bcolors
import colorama
import argparse
from random import random

from norec4dna.Encoder import Encoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.OnlineBPDecoder import OnlineBPDecoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.LTBPDecoder import LTBPDecoder
from norec4dna.LTDecoder import LTDecoder
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.OnlineDecoder import OnlineDecoder
from norec4dna.helper import should_drop_packet

if os.name == "nt" and "PYCHARM_HOSTED" not in os.environ:
    colorama.init()
useDNARules = True
DNARules = FastDNARules()


def blackbox(encoder, decoder):
    encoder.encode_to_packets()
    encoded_packets = encoder.get_encoded_packets()  # list(itertools.permutations(..))
    i = len(encoded_packets)
    invalid_drop = 0
    print(bcolors.BOLD + "[+] Created " + str(i) + " Packets" + bcolors.ENDC)
    for packet in encoded_packets:
        i -= 1
        if isinstance(decoder, LTBPDecoder) or isinstance(decoder, OnlineBPDecoder):
            if not should_drop_packet(DNARules, packet):
                if decoder.input_new_packet(packet):
                    decoder.solve()
                    print(
                        bcolors.OK
                        + "[!] Blackbox dropped "
                        + str(i)
                        + " random Packets, DNA-Simulator dropped "
                        + str(invalid_drop)
                        + " Packets and finished successful"
                        + bcolors.ENDC
                    )
                    return decoder.is_decoded(), len(encoded_packets), i
            else:
                invalid_drop += 1
        else:
            # we are in the FAST_decoder part...
            if not should_drop_packet(DNARules, packet):
                decoder.input_new_packet(packet)
            else:
                invalid_drop += 1
            # only try to solve if we
            if (
                    packet.total_number_of_chunks <= len(encoded_packets) - i - invalid_drop
            ) and decoder.solve():
                print(
                    bcolors.OK
                    + "[!] Blackbox dropped "
                    + str(i)
                    + " random Packets, DNA-Simulator dropped "
                    + str(invalid_drop)
                    + " Packets and finished successful"
                    + bcolors.ENDC
                )
                return decoder.is_decoded(), len(encoded_packets), i
    if decoder.GEPP is not None:
        decoder.solve()

    j = 0  # Make it terminate hard if we doubled our # of encoded Packets...
    while (not (decoder.is_decoded() and decoder.solve())) or j >= 2 * len(
            encoded_packets
    ):
        packet = encoder.create_and_add_new_packet()
        if not should_drop_packet(DNARules, packet):
            decoder.input_new_packet(packet)
        else:
            invalid_drop += 1
        j += 1
    print(
        bcolors.ITALIC
        + "[!] Blackbox dropped "
        + str(i)
        + " random Packets."
        + bcolors.ENDC
    )
    print(
        bcolors.ITALIC
        + "[!] DNA-Simulator dropped "
        + str(invalid_drop)
        + " Packets."
        + bcolors.ENDC
    )
    return decoder.is_decoded(), len(encoder.get_encoded_packets()), i


def __old_should_drop_packet(packet):
    if not useDNARules:
        return False
    rand = random()
    drop_chance = DNARules.apply_all_rules(packet)
    # drop packet if rand bigger than the drop_chance for this Packet.
    return drop_chance > rand


def blackboxOnlineTest(file, number_of_chunks=800, seed=2):
    start = time.time()
    epsilon = 0.05
    quality = 7
    dist = OnlineDistribution(epsilon, seed)
    number_of_chunks = dist.get_size()
    print("Starting Blackbox Test with " + str(number_of_chunks) + " Chunks")
    pseudo = OnlineBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = OnlineEncoder(
        file, number_of_chunks, dist, epsilon, quality, pseudo_decoder=pseudo
    )  # , chunk_size=chunk_size)
    decoder = OnlineDecoder.pseudo_decoder(number_of_chunks)
    decoder.set_read_all_before_decode(True)

    result, numberOfEncodedPackets, droped_count = blackbox(encoder, decoder)
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
        droped_count,
        round(end, 4),
        number_of_chunks,
    ]


def blackboxLTTest(file, number_of_chunks=800, seed=2, chunk_size=0):
    print("Starting Blackbox Test with " + str(number_of_chunks) + " Chunks")
    start = time.time()
    dist = RobustSolitonDistribution(S=number_of_chunks, seed=seed)
    pseudo = LTBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo, chunk_size=chunk_size
    )
    decoder = LTDecoder.pseudo_decoder(number_of_chunks)

    result, numberOfEncodedPackets, dropedCount = blackbox(encoder, decoder)
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
        round(end, 4),
        number_of_chunks,
    ]


def blackboxRU10Test(file, number_of_chunks=800, seed=2, chunk_size=200):
    print("Starting Blackbox Test with " + str(number_of_chunks) + " Chunks")
    start = time.time()
    dist = RaptorDistribution(number_of_chunks)
    pseudo = RU10Decoder.pseudo_decoder()
    encoder = RU10Encoder(
        file, number_of_chunks, dist, pseudo_decoder=pseudo, chunk_size=chunk_size
    )
    decoder = RU10Decoder.pseudo_decoder(number_of_chunks)

    result, numberOfEncodedPackets, droped_count = blackbox(encoder, decoder)
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
        round(end, 4),
        number_of_chunks,
    ]


def get_random_int(max_int):
    return int(random() * max_int)


def main(file="../../main.pdf", mode="LT", repreats=5):
    csv = [
        "filename, codecName, number_of_chunks, numberOfEncodedPackets, droprate, dropped_count, seed, result, time_needed"
    ]
    # "OUT_raptor.pdf" #"b_lq.webm1"  # "b_mq.webm" #
    name = "ERROR"
    # for droprate in np.arange(0.01, 0.02, 0.01):
    for _ in range(repreats):
        rnd = get_random_int(math.pow(2, 31) - 1)
        chunk_size = 20
        insert_header = True
        if chunk_size != 0:
            number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(
                file, chunk_size=chunk_size, insert_header=insert_header
            )
        else:
            number_of_chunks = 800

        try:
            if mode.lower() == "online":
                (
                    name,
                    result,
                    numberOfEncodedPackets,
                    dropped_count,
                    time_needed,
                    number_of_chunks,
                ) = blackboxOnlineTest(
                    file, number_of_chunks=number_of_chunks, seed=rnd
                )
                # chunk_size=chunk_size)
            elif mode.lower() == "lt":
                (
                    name,
                    result,
                    numberOfEncodedPackets,
                    dropped_count,
                    time_needed,
                    number_of_chunks,
                ) = blackboxLTTest(
                    file,
                    number_of_chunks=number_of_chunks,
                    seed=rnd,
                    chunk_size=chunk_size,
                )
            else:
                (
                    name,
                    result,
                    numberOfEncodedPackets,
                    dropped_count,
                    time_needed,
                    number_of_chunks,
                ) = blackboxRU10Test(
                    file,
                    number_of_chunks=number_of_chunks,
                    seed=rnd,
                    chunk_size=chunk_size,
                )

            line = (
                    str(file)
                    + ","
                    + str(name)
                    + ","
                    + str(number_of_chunks)
                    + ","
                    + str(numberOfEncodedPackets)
                    + ","
                    + str(dropped_count)
                    + ","
                    + str(rnd)
                    + ","
                    + str(result)
                    + ","
                    + str(time_needed)
            )
        except Exception as ex:
            print("Error...", ex)
            line = (
                    str(file)
                    + ","
                    + str(name)
                    + ","
                    + str(number_of_chunks)
                    + ","
                    + "ERROR"
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
    dtimeno = "sim" + str(time.strftime("%Y-%m-%d_%H-%M", time.localtime())) + ".csv"
    with open(dtimeno, "w") as f:
        for line in csv:
            f.write(line + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze network traffic")
    parser.add_argument(
        "-f", "--file", help="File to use for Simulation", required=True
    )
    parser.add_argument(
        "-p",
        "--profile",
        help='If set, a profiling Graph "pycallgraph.png" will be created.',
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "-i",
        "--repeats",
        type=int,
        help="Number of repeats the Simulator should run (default=5)",
        default=5,
    )
    parser.add_argument(
        "-m",
        "--mode",
        help="Mode in which the simulator should run. Choose from=[LT(default),Online,Raptor]",
        default="LT",
    )
    parser.add_argument(
        "-d",
        "--dnarules",
        help="If set, DNARules will be applied to all packets.",
        action="store_true",
        required=False,
    )
    args = parser.parse_args()
    profile = bool(args.profile)
    filename = str(args.file)
    mode = str(args.mode)
    useDNARules = bool(args.dnarules)
    repreats = int(args.repeats)
    if profile:
        from pycallgraph import PyCallGraph
        from pycallgraph.output import GraphvizOutput

        print(
            bcolors.WARN
            + "[!] running with profiler - this might decrease performance"
            + bcolors.ENDC
        )
        with PyCallGraph(output=GraphvizOutput()):
            main(filename, mode, repreats)
        print(
            bcolors.BLUE
            + '[*] profiling Graph saved as "pycallgraph.png"'
            + bcolors.ENDC
        )
    else:
        main(filename, mode, repreats)
else:
    print("[!] WARNING: RUN THIS SCRIPT DIRECTLY!")
