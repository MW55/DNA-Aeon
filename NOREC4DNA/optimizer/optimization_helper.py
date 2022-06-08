import copy
import csv
import math
import matplotlib.pyplot as plt
import numpy as np

from norec4dna import Encoder, nocode, RU10Encoder, RU10Decoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper import should_drop_packet
from norec4dna.helper.RU10Helper import intermediate_symbols
from norec4dna.rules.FastDNARules import FastDNARules


def list_to_diff_list(x):
    """
    Converts an accumulated  list to a list of differences.
    :param x: List to convert.
    :return: Converted list.
    """
    XX = np.zeros(len(x) - 1)
    for i in range(1, len(x)):
        XX[i - 1] = x[i] - x[i - 1]
    return XX


def diff_list_to_list(x):
    """
    Converts a list of differences to an accumulated list
    :param x: List to convert.
    :return: Converted list.
    """
    XX = np.zeros(len(x) + 1)
    for i in range(1, len(x) + 1):
        XX[i] = x[i - 1] + XX[i - 1]
    return XX


def scale_to(x, max_num):
    """
    Scales values of a list to the given max number.
    :param x: List to scale.
    :param max_num: Number to scale to.
    :return: Scaled list.
    """
    x = x / (np.max(x) if np.max(x) != 0 else 1.0)
    return x * max_num


def plot_dist(diff_lst):
    """
    Simple plot of a distribution. Use diff_list here.
    :param diff_lst:
    :return:
    """
    plt.plot(diff_lst)
    plt.xlabel("Grade")
    plt.ylabel("Probability")
    plt.show()


def select_dist(dist_list, sel=0.4):
    """
    Selects best distributions from a list that contains tuples (distribution diff_list, avg degree-wise errors,
    summed error) for all distributions of the population.
    :param dist_list: List of distributions.
    :param sel: Percentage of distributions to select, default=0.3.
    :return: List containing the best distributions.
    """
    dist_list.sort(key=lambda x: x[2])
    sel_dist_list = []
    # % parents to select from the population
    pars = math.ceil(len(dist_list) * sel)
    if pars < 2:
        pars = 2
    for x in range(0, pars):
        sel_dist_list.append(dist_list[x])
    return sel_dist_list


def merge_dist_avg(dist_1, dist_2):
    """
    Merge two distributions by adding their degree-wise probabilities and dividing them by 2.
    :param dist_1: Difflist of distribution 1.
    :param dist_2: Difflist of distribution 2.
    :return: Merged difflist.
    """
    x_1 = list_to_diff_list(dist_1[0])
    x_2 = list_to_diff_list(dist_2[0])
    new_dist = [sum(x) / 2 for x in zip(x_1, x_2)]
    return norm_list(new_dist)


def merge_crossover(dist_1, dist_2):
    """
    Takes two distributions and swaps a random number of random grades probabilities.
    :param dist_1: Difflist of distribution 1.
    :param dist_2: Difflist of distribution 2.
    :return: Merged difflist.
    """
    swap_pos = np.random.randint(5, 40, np.random.randint(1, 41))
    for pos in swap_pos:
        tmp = dist_1[pos]
        dist_1[pos] = dist_2[pos]
        dist_2[pos] = tmp
    return [norm_list(dist_1), norm_list(dist_2)]


def merge_best_dist(dist_1, dist_2):
    """
    Merges two distributions by taking the degree-wise probabilities with the lower corresponding error for the new
    distribution. Takes tuples (distribution diff_list, avg degree-wise errors, summed error) for distributions.
    :param dist_1: Tuple of distribution 1.
    :param dist_2: Tuple of distribution 2.
    :return: Merged difflist.
    """
    new_dist = [0] * len(dist_1[0])
    zer_ents = []
    for i in range(0, len(dist_1[0])):
        if dist_1[1][i] < dist_2[1][i]:
            min_prob = max(dist_1[0][i], 0.0)
        else:
            min_prob = max(dist_2[0][i], 0.0)
        if min_prob < 0.0005:
            zer_ents.append(i)
            min_prob = 0
        new_dist[i] = min_prob
    n_zer_cnt = 40 - len(zer_ents)
    # TODO only fill empty probs with average
    if n_zer_cnt < 4:
        ind = np.random.choice(zer_ents, size=4 - n_zer_cnt)
        avg_prob = sum([x if x > 0.0005 else 0 for x in new_dist]) / n_zer_cnt
        for j in ind:
            new_dist[j] = avg_prob
    return norm_list(new_dist)


def mutate_random(x, fac=0.5):
    """
    Multiplies every probability of the distribution with a random number between 1-fac and 1+fac. Use diff_list here.
    :param x: Difflist to mutate.
    :param fac: Maximal mutation factor.
    :return: Mutated difflist.
    """
    mut = np.random.uniform(low=1.0 - fac, high=1.0 + fac, size=40)
    x_new = [x[0] * x[1] for x in zip(x, mut)]
    return norm_list(x_new)


def mutate_deg(dist, fac=0.2, min_grads=5, max_grads=31):
    """
    Mutate between min_degs and max_degs random degrees probabilities with a random factor limited by fac. Use diff_list
    here.
    :param dist: Difflist to mutate.
    :param fac: Maximal mutation factor.
    :param min_grads: Minimal grades to mutate.
    :param max_grads: Maximal grades to mutate.
    :return: Mutated difflist.
    """
    ind_lst = np.random.randint(0, 40, size=np.random.randint(min_grads, max_grads))
    for i in ind_lst:
        dist[i] = dist[i] * np.random.uniform(low=1 - fac, high=1 + fac)
    return dist


def norm_list(lst):
    y = sum(lst)
    return [x / y for x in lst]


def init_population(pop_size, add_raptor=True):
    """
    Initializes the first population with the default raptor distribution and random distributions. Returns a list of
    diff_lists representing the distributions.
    :param pop_size: Number of distributions to generate.
    :param add_raptor: True: Add raptor distribution.
    :return: Population of pop_size distributions.
    """
    pop = []
    if add_raptor:
        # Raptor Distribution
        x = np.asarray(
            [0, 10241, 491582, 712794, 831695, 831695, 831695, 831695, 831695, 831695, 948446, 1032189, 1032189,
             1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189,
             1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189, 1032189,
             1032189, 1032189, 1032189, 1048576])
        x = norm_list(list_to_diff_list(x))
        pop.append(x)
        pop_size -= 1
    # Random Distributions
    for _ in range(0, pop_size):
        x = [1 / 40 for _ in range(0, 40)]
        pop.append(mutate_random(x, 1.0))
    return pop


def compute_population_fitness(pop, fast=False):
    """
    Computes the fitness for each distribution and returns all results.
    :param pop:
    :param fast:
    :return:
    """
    pop_fitness = list()
    for dist in pop:
        avg_err_list, avg_err = compute_cost(dist, diff_list=True, fast=fast)
        pop_fitness.append((dist, avg_err_list[1:], avg_err))
    return pop_fitness


def compute_distribution_fitness(raptor_lst, file_lst, runs=25, chunksize=50):
    """
    Computes the fitness of a distribution with a given fitness function for every given file.
    :param raptor_lst: Raptor list of the distribution to compute the fitness for.
    :param file_lst: List of files to encode with the distribution.
    :param runs: Number of En-/Decoding cycles.
    :param chunksize: Chunksize to use.
    :return:
    """
    overhead_lst = list()
    degree_lst = list()
    for file in file_lst:
        res = encode(file, chunksize, raptor_lst, repeats=runs)
        overhead_lst.append(res[0])
        tmp_degree_lst = np.zeros(41)
        for deg in res[1].keys():
            tmp_degree_lst[deg] = sum(res[1][deg]) / len(res[1][deg])
        degree_lst.append(tmp_degree_lst)
    avg_overhead = sum(overhead_lst) / len(overhead_lst)
    avg_degree_err = []
    for x in range(0, 41):
        tmp_err = 0
        cnt = 0
        for lst in degree_lst:
            if lst[x] != 0:
                tmp_err += lst[x]
                cnt += 1
        if cnt != 0:
            avg_degree_err.append(tmp_err / cnt)
        else:
            avg_degree_err.append(0.0)
    avg_degree_err.pop(0)
    return avg_overhead, avg_degree_err


def encode(file, chunk_size, dist, as_dna=True, repeats=15):
    """
    Encodes the file to packets until the pseudo decoder was able to decode it 'repeats' times with the given chunk size
    and the distribution list.
    :param file: File to encode.
    :param chunk_size: Chunksize to use.
    :param dist: The distribution to calculate the average error and overhead for.
    :param as_dna: If true uses the DNA Rules.
    :param repeats: Number of En-/Decoding cycles.
    :return:
    """
    degree_dict = {}
    overhead_lst = []
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size, insert_header=False)
    distribution = RaptorDistribution(number_of_chunks)
    distribution.f = dist
    distribution.d = [x for x in range(0, 41)]
    if as_dna:
        rules = FastDNARules()
    else:
        rules = None
    encoder = RU10Encoder(file, number_of_chunks, distribution, insert_header=False, rules=rules,
                          error_correction=nocode, id_len_format="H", number_of_chunks_len_format="B",
                          save_number_of_chunks_in_packet=False, mode_1_bmp=False)
    encoder.prepare()
    for _ in range(0, repeats):
        encoder.random_state = np.random.RandomState()
        # print("Master-Seed used: " + str(encoder.random_state.get_state()[1][0]))
        pseudo_decoder = create_pseudo_decoder(encoder.number_of_chunks, distribution)
        needed_packets = 0
        while pseudo_decoder.GEPP is None or not pseudo_decoder.is_decoded():
            needed_packets += 1
            packet = encoder.create_new_packet()
            pseudo_decoder.input_new_packet(packet)
            should_drop_packet(rules, packet)
            if packet.get_degree() not in degree_dict:
                degree_dict[packet.get_degree()] = list()
            degree_dict[packet.get_degree()].append(min(packet.error_prob, 1.0))
        overhead = (needed_packets - encoder.number_of_chunks) / 100.0
        overhead_lst.append(overhead)
    return sum(overhead_lst) / len(overhead_lst), degree_dict


def create_pseudo_decoder(number_of_chunks, distribution):
    """
    Creates a new pseudo decoder for the encoder used to encode the files.
    :param number_of_chunks:
    :param distribution: Distribution used for the encoder.
    :return:
    """
    pseudo_decoder = RU10Decoder.pseudo_decoder(number_of_chunks, False)
    if pseudo_decoder.distribution is None:
        pseudo_decoder.distribution = distribution
        pseudo_decoder.numberOfChunks = number_of_chunks
        _, pseudo_decoder.s, pseudo_decoder.h = intermediate_symbols(number_of_chunks, pseudo_decoder.distribution)
        pseudo_decoder.createAuxBlocks()
    return pseudo_decoder


def compute_cost(dist_lst, c_size_list=None, file_list=None, fast=False, diff_list=False):
    """
    Fitness function that takes a either an accumulated list or a diff_list of a distribution to compute the
    errors for. Diff_lists will be scaled up.
    :param dist_lst: Distribution to compute fitness for.
    :param c_size_list: List of chunk_sizes to compute fitness with.
    :param file_list: List of files to compute fitness with.
    :param fast: Only use on file and chunksize
    :param diff_list: True: Scale up diff_list.
    :return:
    """
    tmp_list = copy.deepcopy(dist_lst)
    if c_size_list is None:
        c_size_list = [50, 75, 100]
        if fast:
            c_size_list = [100]
    if file_list is None:
        file_list = ['Dorn', 'Dorn.tar.gz', 'umr_logo_sw_scaled.png']
        if fast:
            file_list = ['Dorn']
    if diff_list:
        dist_lst = scale_to(diff_list_to_list(dist_lst), 1048576)
    degree_packet_costs = dict()
    n = 0
    for p_tmp in range(45):
        degree_packet_costs[p_tmp] = list()
    for enc_file in file_list:
        for c_size in c_size_list:
            degree_packet_costs1, n1 = encode(enc_file, dist_lst, True, chunk_size=c_size)
            [y.extend(degree_packet_costs1[x]) for x, y in degree_packet_costs.items()]
            n += n1
    n = 1.0 * n / (len(c_size_list) + 1)
    avg_err_per_degree = np.zeros(41)
    for deg in degree_packet_costs.keys():
        if len(degree_packet_costs[deg]) > 0:
            avg_err_per_degree[deg] = sum(degree_packet_costs[deg]) / len(degree_packet_costs[deg])
    avg_err_per_degree[0] -= n
    rel_degs = [0] * len(avg_err_per_degree)
    no_rel_degs = 0
    for x in range(1, len(avg_err_per_degree)):
        if tmp_list[x - 1] >= 0.0001:
            rel_degs[x - 1] = avg_err_per_degree[x - 1]
            no_rel_degs += 1
    summed_error = sum(rel_degs) / no_rel_degs
    return avg_err_per_degree + n, max(0.0, summed_error)


def encode(file, dist_lst, asdna=True, chunk_size=50):
    """

    :param file:
    :param dist_lst:
    :param asdna:
    :param chunk_size:
    :return:
    """
    packets_needed = 0
    packets = dict()
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    dist = RaptorDistribution(number_of_chunks)
    dist.f = dist_lst
    d = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
         29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
    dist.d = d
    dna_rules = FastDNARules()
    if asdna:
        rules = dna_rules
    else:
        rules = None
    x = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=False, rules=rules,
                    error_correction=nocode, id_len_format="H", number_of_chunks_len_format="B",
                    save_number_of_chunks_in_packet=False, mode_1_bmp=False)
    x.prepare()
    y = RU10Decoder.pseudo_decoder(x.number_of_chunks, False)
    if y.distribution is None:  # self.isPseudo and
        y.distribution = RaptorDistribution(x.number_of_chunks)
        y.distribution.f = dist_lst
        y.distribution.d = d
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


def generate_log(file_con, file):
    with open(file + ".csv", 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',', lineterminator='\n', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['gen', 'avg_par_err', 'timestamp', 'dist_lst'])
        for item in file_con:
            writer.writerow([item[0], item[1], item[2], item[3]])


def generate_plot(file_con, file, max_gen):
    avg_errs = [x[1] for x in file_con]
    fig = plt.figure()
    axis1 = fig.add_subplot(211)
    axis1.plot(avg_errs)
    x = [x for x in range(0, max_gen + 1)]
    z = np.polyfit(x, avg_errs, 1)
    p = np.poly1d(z)
    axis1.plot(x, p(x), "r--")
    axis2 = fig.add_subplot(212)
    axis2.plot(file_con[-1][3][0])
    fig.savefig(file + ".png")
    fig.savefig(file + ".svg", format='svg', bbox_inches='tight')


if __name__ == '__main__':
    print()
