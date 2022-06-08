#!/usr/bin/python
# -*- coding: latin-1 -*-
import time
import numpy as np

from norec4dna.LTBPDecoder import LTBPDecoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution


def main():
    file = "../.INFILES/b_lq.webm"

    for number_of_chunks in [500, 600, 700, 800, 900, 1000, 1100, 1200]:
        # for epsilon in [0.03]:
        for overhead in np.arange(0.05, 0.50, 0.01):
            for _ in range(3):
                dist = RobustSolitonDistribution(S=number_of_chunks)
                encoder = LTEncoder(file, number_of_chunks, dist)
                encoder.set_overhead_limit(overhead)
                encoder.encode_to_packets()
                aprrox_decoder = LTBPDecoder(None)
                start = time.time()
                for packet in encoder.encodedPackets:
                    aprrox_decoder.input_new_packet(packet)
                aprrox_decoder.solve()
                end = time.time() - start
                print(
                    "Approx.,"
                    + file
                    + ","
                    + str(number_of_chunks)
                    + ","
                    + str(len(encoder.encodedPackets))
                    + ","
                    + str(aprrox_decoder.is_decoded())
                    + ","
                    + str(end)
                )


if __name__ == "__main__":
    main()
