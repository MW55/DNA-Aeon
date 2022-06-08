import os
import csv
import sys
import time
import signal
import argparse
import matplotlib
import numpy as np
import multiprocessing

from optimization_helper import generate_log, init_population, compute_population_fitness, select_dist, \
    mutate_deg, generate_plot, norm_list, merge_crossover, merge_best_dist

matplotlib.use('Agg')
import matplotlib.pyplot as plt


class EvolutionaryOptimizer:
    def __init__(self, pop_size, max_gen, log, merge, mut_rate, fast=False):
        self.pop_size = pop_size
        self.max_gen = max_gen
        self.log = log
        self.merge = merge
        self.fast = fast
        self.mut_rate = mut_rate
        try:
            if log:
                self.pid = os.getpid()
                self.dir = "EvOpt_" + str(pop_size) + "_" + str(max_gen) + "_" + str(merge) + "_" + str(
                    mut_rate) + "_" + str(self.pid)
                os.mkdir(self.dir)
                self.file = self.dir + "/"
                self.file_con = []
                signal.signal(signal.SIGHUP, self.read_sig)
                signal.signal(signal.SIGKILL, self.kill_sig)
        except FileExistsError:
            print("Dir already exists, using it anyway.")
        except:
            print("Can't use signals.")

    def read_sig(self, signal_number, frame):
        generate_log(self.file_con, "EvOpt_" + str(self.pid) + "cache")
        print("File saved, process continuing.")

    def kill_sig(self, signal_number, frame):
        generate_log(self.file_con, "EvOpt_+" + str(self.pid) + "end")
        print("File saved, process killed.")
        sys.exit()

    def evolutionary_optimization(self, population=None):
        """
        Generate population, compute population fitness, select fittest as parents, merge parent distributions, repeat.
        :param population:
        :return:
        """
        # initialize the first population
        if population is None:
            population = init_population(self.pop_size)
        for gen in range(0, self.max_gen):
            # compute population fitness and select best distributions as parent
            pop_fitness = compute_population_fitness(population, self.fast)
            sel_dist_list = select_dist(pop_fitness)
            population = self.compute_population(sel_dist_list)
            # random mutations
            population = [mutate_deg(x, self.mut_rate, 5, 20) for x in population]
            if self.log:
                avg_par_err = sum(x[2] for x in sel_dist_list) / len(sel_dist_list)
                self.file_con.append([gen, avg_par_err, time.time(), [[x[0], x[2]] for x in sel_dist_list]])
        # get best distribution
        fittest_dist = select_dist(compute_population_fitness(population, self.fast))[0]
        # logging and plotting
        if self.log:
            self.file_con.append([self.max_gen, fittest_dist[2], time.time(), [fittest_dist[0], fittest_dist[2]]])
            generate_log(self.file_con, self.file + "log")
            generate_plot(self.file_con, self.file + "plot", self.max_gen)
        return fittest_dist, self.pid

    def compute_population(self, sel_dist_list):
        """
        Computes a new population from the fittest distributions of the previous one with the choosen merge strategy.
        :param sel_dist_list:
        :return:
        """
        pop = []
        ind_lst = [x for x in range(0, len(sel_dist_list))]
        err_lst = [x[2] for x in sel_dist_list]
        prob_lst = norm_list(err_lst)
        if self.merge == 'crossover':
            dist_lst = [x[0] for x in sel_dist_list]
            for _ in range(0, int(self.pop_size / 2)):
                dist_choice = np.random.choice(ind_lst, size=2, p=prob_lst)
                res = merge_crossover(dist_lst[dist_choice[0]], dist_lst[dist_choice[1]])
                pop.extend(res)
        if self.merge == 'best' or self.merge == 'bestfit':
            if self.merge == 'bestfit':
                pop.append(sel_dist_list[0][0])
            while len(pop) < self.pop_size:
                dist_choice = np.random.choice(ind_lst, size=2, p=prob_lst)
                res = merge_best_dist(sel_dist_list[dist_choice[0]], sel_dist_list[dist_choice[1]])
                pop.append(res)
        if len(pop) > self.pop_size:
            del pop[-1]
        return pop


def main(params, fast=True):
    """
    :param params: Tuple of params: (pop_size, max_gen, log, merge)
    :param fast:
    :return:
    """
    opt = EvolutionaryOptimizer(params[0], params[1], params[2], params[3], params[4], fast)
    res = opt.evolutionary_optimization()
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--pop_size", metavar="pop_size", required=False, type=int, action="append",
                        help="list of population sizes for each generation")
    parser.add_argument("--max_gen", metavar="max_gen", required=False, type=int, action="append",
                        help="list of numbers of generations to compute")
    parser.add_argument("--mut_rate", metavar="mut_rate", required=False, type=float, action="append",
                        help="mutation rate for random mutations")
    parser.add_argument("--log", metavar="log", required=False, type=bool, default=False,
                        help="log results as csv")
    parser.add_argument("--spare1core", required=False, default=False, action="store_true")
    parser.add_argument("--merge", metavar="merge", required=False, type=str, action="append",
                        help='choose between crossover, best and bestfit')
    parser.add_argument("--cores", required=False, type=int)
    args = parser.parse_args()
    filename = args.filename
    pop_size = args.pop_size
    max_gen = args.max_gen
    mut_rate = args.mut_rate
    log = args.log
    spare1core = args.spare1core
    merge = args.merge
    if args.cores:
        cores = args.cores
    else:
        cores = multiprocessing.cpu_count()
    params = []
    if spare1core:
        cores = cores - 1
    # fill missing params
    if len(pop_size) < cores:
        [pop_size.append(2) for _ in range(0, cores - len(pop_size))]
    if len(max_gen) < cores:
        [max_gen.append(2) for _ in range(0, cores - len(max_gen))]
    if len(merge) < cores:
        [merge.append('bestfit') for _ in range(0, cores - len(merge))]
    if len(mut_rate) < cores:
        [mut_rate.append(0.01) for _ in range(0, cores - len(mut_rate))]
    for i in range(0, cores):
        params.append((pop_size[i], max_gen[i], log, merge[i], mut_rate[i]))
    p = multiprocessing.Pool(cores)
    a = p.map(main, params)
    # log and plot overall results
    if log:
        with open("EvOpt_overall.csv", 'w', newline='') as f:
            writer = csv.writer(f, delimiter=',', lineterminator='\n', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['dist', 'errs', 'avg_err'])
            for item in a:
                writer.writerow([item[0][0], item[0][1], item[0][2]])
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])
        errs = [x[0][2] for x in a]
        labs = [x[1] for x in a]
        ax.bar(labs, errs)
        plt.savefig("EvOpt_overall.png")
