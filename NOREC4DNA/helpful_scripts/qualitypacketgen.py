import matplotlib.pyplot as plt
from numpy import mean

from norec4dna import RU10Encoder, reed_solomon_encode, Encoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper import should_drop_packet
from norec4dna.rules.FastDNARules import FastDNARules


class QualityPacketGen:
    def __init__(self, encoder):
        encoder.prepare()
        i = 0
        tmp_list = []

        while i < 1000:
            packet = encoder.create_new_packet()
            p_res = should_drop_packet(rules, packet)
            tmp_list.append(packet)
            i += 1

        QualityAnalyzer(tmp_list)


class QualityAnalyzer:
    def __init__(self, tmp_list):
        tmp_dict = dict()
        for packet in tmp_list:
            err_prob = packet.error_prob
            for chunk_no in packet.get_used_packets():
                if chunk_no not in tmp_dict:
                    tmp_dict[chunk_no] = []
                tmp_dict[chunk_no].append(err_prob)

        # for key, data_list in tmp_dict.items():
        num_list = []
        mean_list = []
        for key in tmp_dict:
            mean_list.append(len(tmp_dict[key]) / mean(tmp_dict[key]))
            num_list.append(len(tmp_dict[key]))

        index = []
        data = []
        for i, (key, val) in enumerate(sorted(tmp_dict.items())):
            index.append(key)
            data.append(val)

        fig, (ax) = plt.subplots(ncols=1)
        ax.boxplot(data)
        ax.set_xticklabels(index)
        plt.show()

        plt.plot(num_list)
        plt.plot(mean_list)
        plt.show()


if __name__ == "__main__":
    file = "../.INFILES/Dorn"
    chunk_size = 100
    norepairsymbols = 6
    save_number_of_chunks_in_packet = False
    insert_header = False
    rules = FastDNARules()
    error_correction = lambda x: reed_solomon_encode(x, norepairsymbols)
    number_of_chunks = 50
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)

    dist = RaptorDistribution(number_of_chunks)
    x = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=insert_header, rules=rules,
                    error_correction=error_correction, id_len_format="H", number_of_chunks_len_format="B",
                    save_number_of_chunks_in_packet=save_number_of_chunks_in_packet)
    aa = QualityPacketGen(x)
