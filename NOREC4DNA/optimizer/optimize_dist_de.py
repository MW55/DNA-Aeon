import os
import time
import argparse
import numpy as np
import multiprocessing

from optimization_helper import init_population, compute_population_fitness, select_dist, generate_log, \
    generate_plot, compute_cost

"""if self.merge == 'diff':
    
    """


class DifferentialOptimizer:
    def __init__(self, pop_size, max_gen, cr, f, log, fast):
        self.pop_size = pop_size
        self.max_gen = max_gen
        self.cr = cr
        self.f = f
        self.log = log
        self.fast = fast
        if log:
            self.pid = os.getpid()
            self.dir = "DEOpt_" + str(pop_size) + "_" + str(max_gen) + "_" + str(cr) + "_" + str(f) + "_" + str(
                self.pid)
            os.mkdir(self.dir)
            self.file = self.dir + "/"
            self.file_con = []

    def differential_optimization(self, population=None):
        if population is None:
            population = init_population(self.pop_size)
        for gen in range(0, self.max_gen):
            # compute population fitness and new population
            pop_fitness = compute_population_fitness(population, self.fast)
            population = self.compute_population_diff(pop_fitness)
            if self.log:
                avg_par_err = sum(x[2] for x in population) / len(population)
                self.file_con.append([gen, avg_par_err, time.time(), [[x[0], x[2]] for x in population]])
            fittest_dist = select_dist(population)[0]
            population = [x[0] for x in population]
        if self.log:
            self.file_con.append([self.max_gen, fittest_dist[2], time.time(), [fittest_dist[0], fittest_dist[2]]])
            generate_log(self.file_con, self.file + "log")
            generate_plot(self.file_con, self.file + "plot", self.max_gen)
        return fittest_dist

    def compute_population_diff(self, pop_fitness):
        pop = []
        for dist in pop_fitness:
            new_dist = dist[0]
            tmp_dist = [d for d in pop_fitness if not np.array_equal(d, dist)]
            ind = np.random.randint(0, len(tmp_dist), size=3)
            rand = []
            for x in ind:
                rand.append(tmp_dist[x])
            probs = np.random.uniform(0.0, 1.0, size=40)
            for i in range(0, len(probs)):
                if probs[i] < self.cr:
                    new_dist[i] = rand[0][0][i] + self.f * (rand[1][0][i] - rand[2][0][i])
                    if new_dist[i] < 0:
                        new_dist[i] = 0
            new_dist = new_dist / sum(new_dist)
            err_list, err_sum = compute_cost(new_dist, diff_list=True, fast=self.fast)
            new_dist = (new_dist, err_list, err_sum)
            if new_dist[2] < dist[2]:
                pop.append(new_dist)
            else:
                pop.append(dist)
        return pop


def main(params, fast=True):
    """
    :param params: Tuple of params: (pop_size, max_gen, log, merge)
    :param fast:
    :return:
    """
    opt = DifferentialOptimizer(params[0], params[1], params[2], params[3], params[4], fast)
    res = opt.differential_optimization()
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--pop_size", metavar="pop_size", required=False, type=int, action="append",
                        help="list of population sizes for each generation")
    parser.add_argument("--max_gen", metavar="max_gen", required=False, type=int, action="append",
                        help="list of numbers of generations to compute")
    parser.add_argument("--cr", metavar="mut_rate", required=False, type=float, action="append",
                        help="crossover rate")
    parser.add_argument("--f", metavar="f", required=False, type=float, action="append",
                        help="differential weight")
    parser.add_argument("--log", metavar="log", required=False, type=bool, default=False,
                        help="log results as csv")
    parser.add_argument("--spare1core", required=False, default=False, action="store_true")
    parser.add_argument("--cores", required=False, type=int)
    args = parser.parse_args()
    filename = args.filename
    pop_size = args.pop_size
    max_gen = args.max_gen
    cr = args.cr
    f = args.f
    log = args.log
    spare1core = args.spare1core
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
    if len(cr) < cores:
        [cr.append(0.3) for _ in range(0, cores - len(cr))]
    if len(f) < cores:
        [f.append(0.5) for _ in range(0, cores - len(f))]
    for i in range(0, cores):
        params.append((pop_size[i], max_gen[i], cr[i], f[i], log))
    p = multiprocessing.Pool(cores)
    a = p.map(main, params)
    # main((4, 5, 0.4, 0.5, True), True)
