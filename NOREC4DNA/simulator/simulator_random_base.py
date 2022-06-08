import random
import argparse
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

import norec4dna.rules.DNARules


class RandomBaseSimulator:
    def __init__(self, base_length):
        self.length = base_length

    def generate_random_base_string(self):
        res = ""
        for _ in range(self.length):
            res += random.choice(["A", "C", "G", "T"])
        return res

    def run_simulation(self, repeats=5, out_file=None):
        rule_obj = norec4dna.rules.FastDNARules.FastDNARules()
        res = []

        for _ in range(repeats):
            drop_chance, data, pkt = rule_obj.apply_all_rules_with_data(self.generate_random_base_string())
            res.append(drop_chance)
        if out_file is not None:
            with open(out_file, "w") as f:
                for csv_line in res:
                    f.write(csv_line + "\n")
        return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze DNA-Rules with random Basestrings"
    )
    parser.add_argument(
        "-o", "--outfile", help="Name of the Outfile", required=False, default=None
    )
    parser.add_argument(
        "-l",
        "--baselength",
        type=int,
        help="Length of each BaseString (default=400)",
        default=200,
    )
    parser.add_argument(
        "-i",
        "--repeats",
        type=int,
        help="Number of repeats the Simulator should run (default=5)",
        default=5,
    )
    args = parser.parse_args()
    filename = args.outfile
    repeats = int(args.repeats)
    length = int(args.baselength)

    simulator = RandomBaseSimulator(length)
    err_prob_list = simulator.run_simulation(repeats, filename)
    binwidth = 0.01
    hist_dist = scipy.stats.rv_histogram(np.histogram(err_prob_list, bins=20))
    X = np.linspace(0, 4, 50)
    plt.xlabel("Error probability as calculated by the ruleset")  # can be >1.0 because it is additive!
    plt.ylabel("Density")
    plt.plot(X, hist_dist.pdf(X))
    plt.plot(X, hist_dist.cdf(X))
    plt.grid(True)
    plt.savefig("hist_pdf_cdf.pdf")
    plt.show()
