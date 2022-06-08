import numpy as np
import time, math
from random import random
from datetime import datetime

from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.OnlineBPDecoder import OnlineBPDecoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.LTBPDecoder import LTBPDecoder
from norec4dna.LTDecoder import LTDecoder
from norec4dna.distributions.IdealSolitonDistribution import IdealSolitonDistribution
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution
from norec4dna.distributions.OnlineDistribution import OnlineDistribution


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
    print(
        "Finished encoding after "
        + str(round(end, 4))
        + " sec. "
        + str(len(encoder.get_encoded_packets()))
        + " Packets encoded.\n"
    )

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


""" encoder -> blackbox with a (droprate) change of dropping packet -> decoder """


def blackbox(encoder, decoder, droprate=0.02):
    encoder.encode_to_packets()
    i = 0
    print("[+] Created " + str(len(encoder.get_encoded_packets())) + " Packets")
    for packet in encoder.get_encoded_packets():
        if np.random.rand() > droprate:
            decoder.input_new_packet(packet)
        else:
            i += 1
    print("[!] Blackbox dropped " + str(i) + " random Packets.")
    return (decoder.is_decoded(), len(encoder.get_encoded_packets()), i)


def blackboxTest(file, number_of_chunks=800, droprate=0.02, seed=2):
    print(
        "#### Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks and a droprate of "
        + str(droprate)
        + " ####"
    )
    start = time.time()
    dist = IdealSolitonDistribution(number_of_chunks, seed)
    pseudo = LTBPDecoder.pseudo_decoder(number_of_chunks)
    encoder = LTEncoder(file, number_of_chunks, dist, pseudo_decoder=pseudo)
    decoder = LTDecoder.pseudo_decoder(number_of_chunks)

    result, numberOfEncodedPackets, dropedCount = blackbox(
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
    return ["LT", result, numberOfEncodedPackets, dropedCount, round(end, 4)]


def get_random_int(maxInt):
    return int(random() * maxInt)


if __name__ == "__main__":
    dtimenow = datetime.now().strftime("%Y-%m-%d_%H-%M")
    csv = [
        "filename, codecName, number_of_chunks, numberOfEncodedPackets, droprate, dropedCount, seed, result, timeNeeded"
    ]
    file = "b_lq.webm1"  # "b_mq.webm"
    # main(file)
    name = "ERROR"
    for droprate in np.arange(0.01, 0.06, 0.01):
        for repeat in range(10):
            try:
                print("GOGOGO...")
                rnd = get_random_int(math.pow(2, 31) - 1)
                _number_of_chunks = 800
                print("GOGOGO...")
                name, result, numberOfEncodedPackets, dropedCount, timeNeeded = blackboxTest(
                    file, number_of_chunks=_number_of_chunks, droprate=droprate, seed=rnd
                )
                line = (
                        str(file)
                        + ", "
                        + str(name)
                        + ", "
                        + str(_number_of_chunks)
                        + ", "
                        + str(numberOfEncodedPackets)
                        + ", "
                        + str(droprate)
                        + ", "
                        + str(dropedCount)
                        + ", "
                        + str(rnd)
                        + ", "
                        + str(result)
                        + ", "
                        + str(timeNeeded)
                )
            except (Exception):
                print("Error...")
                line = (
                        str(file)
                        + ", "
                        + str(name)
                        + ", "
                        + str(_number_of_chunks)
                        + ", "
                        + "ERROR"
                        + ", "
                        + str(droprate)
                        + ", "
                        + "ERROR"
                        + ", "
                        + str(rnd)
                        + ", "
                        + "ERROR"
                        + ", "
                        + "ERROR"
                )
            print(line)
            csv.append(line)
    dtimenow = datetime.now().strftime("%Y-%m-%d_%H-%M")
    with open("sim" + str(dtimenow) + ".csv", "w") as f:
        for line in csv:
            f.write(line + "\n")
else:
    print("[!] WARNING: RUN THIS SCRIPT DIRECTLY!")
