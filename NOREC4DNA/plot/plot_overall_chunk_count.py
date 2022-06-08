import typing
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

from norec4dna.Packet import Packet
from norec4dna.RU10Packet import RU10Packet
from norec4dna.OnlinePacket import OnlinePacket
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna import RU10Decoder, OnlineDecoder, RU10Encoder
from norec4dna.ErrorCorrection import get_error_correction_encode
from norec4dna.distributions.RaptorDistribution import RaptorDistribution


def __handle_raptor(packet: RU10Packet):
    p_decoder = RU10Decoder.pseudo_decoder(packet.get_total_number_of_chunks())
    p_decoder.input_new_packet(packet)
    return p_decoder.removeAndXorAuxPackets(packet)


def __handle_online(packet: OnlinePacket):
    p_decoder = OnlineDecoder.pseudo_decoder(packet.get_total_number_of_chunks())
    p_decoder.input_new_packet(packet)
    return p_decoder.removeAndXorAuxPackets(packet)


def __handle_lt(packet: Packet):
    return packet.get_used_packets()


def extract_used(packet: Packet):
    if isinstance(packet, RU10Packet):
        return __handle_raptor(packet)
    elif isinstance(packet, OnlinePacket):
        return __handle_online(packet)
    return __handle_lt(packet)


def plot_chunks_count_in_packets(packets: typing.List[Packet], save_prefix="logo"):
    counter = [0] * packets[0].get_total_number_of_chunks()
    for packet in packets:
        i = 0
        cleared = extract_used(packet)
        for i in range(len(cleared)):
            counter[i] += cleared[i]
    save_str = save_prefix + "_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    with open(save_str + ".csv", "w") as out_file:
        out_file.write(",".join([str(x) for x in counter]))
    plot_numbers(counter, save_str)


def plot_from_csv(file_name: str):
    """
    We cloud just use pandas to plot the csv but this approach splits reading and plotting...
    """
    csv = pd.read_csv(file_name, header=None)
    lst = csv.values.tolist()[0]
    print(lst)
    plot_numbers(lst, file_name.rsplit(".", 1)[0])


def plot_numbers(counter: typing.List[int], save_str):
    plt.plot(counter)
    plt.xlabel("Chunk")
    plt.ylabel("# Occurrences")
    plt.savefig(save_str + '.svg')


if __name__ == "__main__":
    file = "logo.jpg"
    number_of_chunks = 400
    insert_header = True
    rules = FastDNARules()
    norepair_symbols = 4
    error_correction_str = "reedsolomon"
    error_correction = get_error_correction_encode(error_correction_str, norepair_symbols)
    PACKET_LEN_FORMAT = "I"
    CRC_LEN_FORMAT = "I"
    NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
    ID_LEN_FORMAT = "I"
    save_number_of_chunks_in_packet = False
    upper_bound = 0.9
    encoder = RU10Encoder(file, number_of_chunks, RaptorDistribution(number_of_chunks=number_of_chunks),
                          insert_header=insert_header, pseudo_decoder=None,
                          chunk_size=0, rules=rules, error_correction=error_correction,
                          packet_len_format=PACKET_LEN_FORMAT,
                          crc_len_format=CRC_LEN_FORMAT, number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
                          id_len_format=ID_LEN_FORMAT, save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
                          mode_1_bmp=False, prepend="", append="", drop_upper_bound=upper_bound)
    encoder.set_overhead_limit(0.40)
    # encoder.encode_to_packets()
    save_prefix = file + "_" + str(number_of_chunks) + "_" + ("H" if insert_header else "") + "_" + \
                  ("R" if rules is not None else "") + "_" + error_correction_str + "_" + str(norepair_symbols) + \
                  "_" + PACKET_LEN_FORMAT + "_" + CRC_LEN_FORMAT + "_" + NUMBER_OF_CHUNKS_LEN_FORMAT + "_" + \
                  ID_LEN_FORMAT + "_" + ("NUMC" if save_number_of_chunks_in_packet else "") + "_" + str(upper_bound)
    # plot_chunks_count_in_packets(list(encoder.encodedPackets), save_prefix)
    plot_from_csv("logo.jpg_400_H_R_reedsolomon_4_I_I_I_I__0.9_2021-01-20_12-40-52.csv")
