import string
import random
import pickle
import argparse
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from norec4dna.helper.RU10Helper import intermediate_symbols
from norec4dna import Encoder, nocode, RU10Encoder, RU10Decoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper import should_drop_packet
from optimization_helper import list_to_diff_list, scale_to, diff_list_to_list

DO_PLOT = False
DO_MOVIE = True
FPS = 30


class GradientDescentOptimizer:
    def __init__(self, runs, x, d, f_name='Dorn'):
        self.runs = runs
        self.X = x
        self.d = d
        self.x_sum = x[-1]  # sum(X)
        self.FILENAME = f_name
        self.current_min = 100000.42

    def encode(self, file, asdna=True, error_correction=nocode, insert_header=False,
               save_number_of_chunks_in_packet=False, mode_1_bmp=False, chunk_size=50):
        packets_needed = 0
        packets = dict()
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
        dist = RaptorDistribution(number_of_chunks)
        dist.f = self.X
        dist.d = self.d
        dna_rules = FastDNARules()
        if asdna:
            rules = dna_rules
        else:
            rules = None
        x = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=insert_header, rules=rules,
                        error_correction=error_correction, id_len_format="H", number_of_chunks_len_format="B",
                        save_number_of_chunks_in_packet=save_number_of_chunks_in_packet, mode_1_bmp=mode_1_bmp)
        x.prepare()
        y = RU10Decoder.pseudo_decoder(x.number_of_chunks, False)
        if y.distribution is None:  # self.isPseudo and
            y.distribution = RaptorDistribution(x.number_of_chunks)
            y.distribution.f = self.X
            y.distribution.d = self.d
            y.number_of_chunks = x.number_of_chunks
            _, y.s, y.h = intermediate_symbols(x.number_of_chunks, y.distribution)
            y.createAuxBlocks()
        n = 0
        for p_tmp in range(45):
            packets[p_tmp] = list()
        while n < number_of_chunks * 50:
            pack = x.create_new_packet()
            if packets_needed == 0:
                y.input_new_packet(pack)
            should_drop_packet(dna_rules, pack)

            if pack.get_degree() not in packets:
                packets[pack.get_degree()] = list()
            packets[pack.get_degree()].append(pack.error_prob)
            n += 1
            if n >= number_of_chunks and y.is_decoded() and packets_needed == 0:
                packets_needed = n
                # we dont want to break, we want to generate #chunks * XXX packets!
                # break
        print("Packets created: " + str(sum([len(x) for x in packets.values()])))
        return packets, (packets_needed - number_of_chunks) / 100.0

    def compute_cost(self, c_size_list=None, file_list=None):
        if c_size_list is None:
            c_size_list = [50, 75, 100]
        if file_list is None:
            file_list = ['Dorn', 'Dorn.tar.gz', 'umr_logo_sw_scaled.png']
        degree_packet_costs = dict()
        n = 0
        for p_tmp in range(45):
            degree_packet_costs[p_tmp] = list()
        for enc_file in file_list:
            for c_size in c_size_list:
                degree_packet_costs1, n1 = self.encode(enc_file, True, nocode, False, False, False, chunk_size=c_size)
                [y.extend(degree_packet_costs1[x]) for x, y in degree_packet_costs.items()]
                n += n1
        n = 1.0 * n / (len(c_size_list) + 1)
        tmp_list = np.array([len(x) for x in degree_packet_costs.values()])
        avg_err_per_degree = np.zeros(max(self.d) + 1)
        used_degrees = 0
        for deg in degree_packet_costs.keys():
            if len(degree_packet_costs[deg]) > 0:
                avg_err_per_degree[deg] = sum(degree_packet_costs[deg]) / len(degree_packet_costs[deg])
                used_degrees += 1
        avg_err_per_degree[0] -= n
        summed_error = sum(avg_err_per_degree + n) / used_degrees
        if summed_error <= self.current_min:
            if DO_PLOT:
                plt.plot(tmp_list / sum(tmp_list))
                plt.title("Dotted: Distribution - Line: Observed Dist - Err:" + str(summed_error))
            # plt.show(block=True)
            tmp_list = np.insert(list_to_diff_list(self.X), 0, 0)
            print(
                "Created packet (" + str(summed_error) + ") is BETTER than old minimum (" + str(self.current_min) + ")")
            print("Distribution:")
            print(self.X)
            if DO_PLOT:
                plt.plot(tmp_list / sum(tmp_list), ":")
                plt.xticks(np.arange(0, 45, 2))
                plt.grid()
                plt.show(block=False)
            self.current_min = summed_error
        else:
            print(
                "Created packet (" + str(summed_error) + ") is worse than old minimum (" + str(self.current_min) + ")")
        return avg_err_per_degree + n - summed_error, summed_error

    def gradient_descent(self, X, y, alpha):
        cost = []  # np.zeros(self.runs)
        solutions = []
        X = list_to_diff_list(X)
        for i in range(self.runs):
            # h = X @ theta.T
            # loss = h - y
            # gradient = np.dot(X.T, loss) / (i + 1)  # self.runs
            # theta = theta - alpha * gradient
            # theta = theta - (alpha / len(X)) * np.sum(X * (X @ theta.T - y), axis=0)
            theta = alpha * (self.runs - i) * (
                    y - (1.0 - (self.runs - i) / self.runs))  # decrease the stepsize with each run
            X = np.absolute(X - theta)  # / sum(self.X)
            X = np.array(X) / sum(X)
            self.X = np.round(scale_to(diff_list_to_list(X), self.x_sum))
            y, summed_error = self.compute_cost()
            y = y[1:]
            solutions.append(self.X)
            cost.append((y, summed_error))
            print("Summed cost: " + str(np.sum(y)) + " , averaged: " + str(summed_error))
        return solutions, cost

    # running the gd and cost function
    # g, cost = gradientDescent(X, y, theta, iters, alpha)


def main(gd_runs=100):
    # if __name__ == "__main__":
    # f = [0, 10241, 491582, 712794, 831695, 948446, 1032189, 1048576]
    d = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
         30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
    X = [0, 10241, 491582, 712794, 831695, 831695, 831695, 831695, 831695, 831695, 948446, 1032189, 1032189, 1032189,
         1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189,
         1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189,
         1032189, 1032189, 1048576]
    opt = GradientDescentOptimizer(gd_runs, X, d)  # , f_name)
    X = np.array(X) / sum(X)
    """
    # Get current size
    fig_size = plt.rcParams["figure.figsize"]
    fig_dpi = plt.rcParams["figure.dpi"]
    plt.rcParams["svg.fonttype"] = "none"
    # Prints: [8.0, 6.0]
    print("Current size:", fig_size)
    print("Current dpi:", fig_dpi)
    # Set figure width to 12 and height to 9
    fig_size[0] = 16
    fig_size[1] = 9
    fig_dpi = 500
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams["figure.dpi"] = fig_dpi
    """
    res = []
    for _ in range(1):
        y = np.ones(len(X) - 1)
        alpha = np.random.uniform(0.0001, 0.25, 1)
        new_x, tmp = opt.gradient_descent(X, y, alpha)
        for x in tmp:
            if DO_PLOT:
                plt.plot(x[0])
        min_pos = np.argmin([y for x, y in tmp])
        print("Min: ", sum(tmp[min_pos]))
        print("Distr: ")
        print(new_x[min_pos])
        if DO_PLOT:
            plt.show(block=False)
            plt.plot([sum(i) for i, j in tmp])
            plt.title("Error per GD-RUN:")
        res.append((new_x, tmp))
        if DO_MOVIE:
            writer = animation.FFMpegWriter(fps=2 * FPS, metadata=dict(artist='Michael Schwarz'), bitrate=3600)
            fig = plt.figure()
            l, = plt.plot([], [])  # , 'k-o')

            plt.xlim(0, len(X) + 1)
            plt.ylim(0, 1.5)
            with open('reses', 'wb') as picke_out:
                pickle.dump(res, picke_out)
            rand_str = ''.join(random.choice(string.ascii_lowercase) for _ in range(12))
            with writer.saving(fig, "tmp/dist_development_gd_" + rand_str + ".mp4", 500):
                for xx in res:
                    for x in xx[0]:
                        x = list_to_diff_list(scale_to(x, 1.0))
                        l.set_data([i for i in range(len(x))], x)
                        for y in range(FPS):
                            writer.grab_frame()
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    # parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
    #                    help="Error Correction Method to use; possible values: \
    #                    nocode, crc, reedsolomon (default=nocode)")
    parser.add_argument("--runs", metavar="runs", required=False, type=int, default=100,
                        help="number of runs for GradientDecent")
    parser.add_argument("--spare1core", required=False, default=False, action="store_true")
    parser.add_argument("--plot", required=False, default=False, action="store_true")
    args = parser.parse_args()

    filename = args.filename
    runs = args.runs
    spare1core = args.spare1core
    DO_PLOT = args.plot
    cores = multiprocessing.cpu_count()
    if spare1core:
        cores = cores - 1
    p = multiprocessing.Pool(cores)
    param = [None] * cores
    a = p.map(main, [runs for _ in range(cores)])
    X = list(y for x in a for y in x[0][0])
    tmp = list(y for x in a for y in x[0][1])
    print("Absolute min:")
    min_pos = np.argmin([y for x, y in tmp])
    print("Min: ", tmp[min_pos][1])
    print("Distr: ")
    print(X[min_pos])
    if DO_PLOT:
        plt.plot(X[min_pos])
        plt.title("x with min:")
        plt.show(block=False)
        plt.plot([0.] + list_to_diff_list(tmp[min_pos][0]))
        plt.title("tmp:")
        plt.show(block=True)
    with open('tmp/aes', 'wb') as picke_out:
        pickle.dump(a, picke_out)
    with open('tmp/Xes', 'wb') as picke_out:
        pickle.dump(X, picke_out)
    with open('tmp/tmpes', 'wb') as picke_out:
        pickle.dump(tmp, picke_out)

# current config: 100 * 3 * 3 * #cores
# server: 72000 encoder runs with ~[1550 2250 2500] packets each
# yields: [111600000, 162000000, 180000000] packets = 453600000 generated overall
