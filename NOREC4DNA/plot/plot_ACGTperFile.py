import typing
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

font = {"family": "normal", "size": 15}
axis_font = {"size": "20"}
matplotlib.rc("font", **font)
plt.rcParams["svg.fonttype"] = "none"


# Ähnlich wie in "A highly parallel strategy for storage of
# digital information in living cells"
# https://www.biorxiv.org/content/biorxiv/early/2016/12/26/096792.full.pdf


def main(text, out_filename, group: int = 32):
    sumA = None
    sumT = None
    sumG = None
    sumC = None
    fi = 1
    j = 0
    print(text)
    leng = np.int(np.ceil(len(text) / group))
    A = np.zeros(leng, dtype=np.float64)
    C = np.zeros(leng, dtype=np.float64)
    G = np.zeros(leng, dtype=np.float64)
    T = np.zeros(leng, dtype=np.float64)
    error = 0
    for i in range(0, len(text), group):
        tx = text[i: i + group]
        for t in tx:
            if t == "A":
                A[j] += 1
            elif t == "C":
                C[j] += 1
            elif t == "G":
                G[j] += 1
            elif t == "T":
                T[j] += 1
            else:
                error += 1
        A[j] = 1.0 * A[j] / (1.0 * len(tx))
        C[j] = 1.0 * C[j] / (1.0 * len(tx))
        G[j] = 1.0 * G[j] / (1.0 * len(tx))
        T[j] = 1.0 * T[j] / (1.0 * len(tx))

        j += 1
    print(A)
    plt.plot(A, label="A")
    plt.plot(C, label="C")
    plt.plot(G, label="G")
    plt.plot(T, label="T")
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
    plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(xmax=j - 1))
    plt.ylabel("relative density")
    plt.xlabel("generated sequence")
    plt.grid(True)
    plt.legend()
    plt.show(block=False)
    plt.savefig(out_filename + "_" + str(group) + ".pdf", bbox_inches="tight")
    plt.savefig(out_filename + "_" + str(group) + ".svg")
    plt.close()
    # Boxplot:
    plt.boxplot([A, C, G, T])  # , label="A")
    plt.xticks([1, 2, 3, 4], ["A", "C", "G", "T"])
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    plt.tight_layout()
    plt.grid(True)
    plt.show(block=False)
    plt.savefig(out_filename + "_" + "box_" + str(group) + ".pdf", bbox_inches="tight")
    plt.savefig(out_filename + "_" + "box_" + str(group) + ".svg")
    plt.close()

    if sumA is None:
        sumA = A
        sumC = C
        sumG = G
        sumT = T
    else:
        sumA = np.add(sumA, A)
        sumC = np.add(sumC, C)
        sumG = np.add(sumG, G)
        sumT = np.add(sumT, T)

    if sumA is not None:
        sumA = np.divide(sumA, 1.0 * fi)
        sumC = np.divide(sumC, 1.0 * fi)
        sumG = np.divide(sumG, 1.0 * fi)
        sumT = np.divide(sumT, 1.0 * fi)
        print(sumA)
        plt.plot(sumA, label="A")
        plt.plot(sumC, label="C")
        plt.plot(sumG, label="G")
        plt.plot(sumT, label="T")
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
        plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(xmax=j - 1))
        plt.ylabel("relative density")
        plt.xlabel("generated sequence")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.show(block=False)
        print(out_filename + "_average_" + str(group) + ".pdf")
        plt.savefig(out_filename + "_average_" + str(group) + ".pdf", bbox_inches="tight", )
        plt.savefig(out_filename + "_average_" + str(group) + ".svg")
        plt.close()


if __name__ == "__main__":
    filename = ".OUTFILES/output.fasta"
    with open(filename, "r") as in_file:
        lines: typing.List[str] = in_file.readlines()
    text = lines[1]
    main(text, filename)
