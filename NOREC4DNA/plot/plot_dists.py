import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("..")
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution
from norec4dna.distributions.IdealSolitonDistribution import IdealSolitonDistribution
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution

plt.rcParams["svg.fonttype"] = "none"


def vergleich():
    S = 50
    for a in ["log", ""]:
        delt = 0.1
        robust = RobustSolitonDistribution(S=S, K=8, delta=delt, seed=0)
        print(robust.pre_comp_dist)
        if a == "log":
            plt.semilogy([0.0] + robust.pre_comp_dist, label="K=8,   $\delta$ = 0.1")
        else:
            plt.plot([0.0] + robust.pre_comp_dist, label="K=8,   $\delta$ = 0.1")

        robust = RobustSolitonDistribution(S=S, K=8, delta=1.0, seed=0)
        print(robust.pre_comp_dist)
        if a == "log":
            plt.semilogy([0.0] + robust.pre_comp_dist, label="K=8,   $\delta$ = 1.0")
        else:
            plt.plot([0.0] + robust.pre_comp_dist, label="K=8,   $\delta$ = 1.0")

        robust = RobustSolitonDistribution(S=S, K=15, delta=0.5, seed=0)
        print(robust.pre_comp_dist)
        if a == "log":
            plt.semilogy([0.0] + robust.pre_comp_dist, label="K=15, $\delta$ = 0.5")
        else:
            plt.plot([0.0] + robust.pre_comp_dist, label="K=15, $\delta$ = 0.5")

        plt.ylabel("Probability")
        plt.xlabel("Degree")
        manager = plt.get_current_fig_manager()
        # manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.show(block=False)
        plt.savefig("../plotDists/Vergleich_RobustSolitonS50_K8_delta" + str(delt) + "_" + a + ".pdf",
                    bbox_inches="tight", )
        plt.savefig("../plotDists/Vergleich_RobustSolitonS50_K8_delta" + str(delt) + "_" + a + ".svg")
        plt.close()


def main():
    S = 50
    for a in ["log", ""]:
        for delt in np.arange(0.1, 1.1, 0.1):
            robust = RobustSolitonDistribution(S=S, K=8, delta=delt, seed=0)
            print(robust.pre_comp_dist)
            if a == "log":
                plt.semilogy([0.0] + robust.pre_comp_dist)
            else:
                plt.plot([0.0] + robust.pre_comp_dist)
            plt.ylabel("Probability")
            plt.xlabel("Degree")
            manager = plt.get_current_fig_manager()
            # manager.resize(*manager.window.maxsize())
            plt.grid(True)
            plt.tight_layout()
            plt.show(block=False)
            plt.savefig("../plotDists/RobustSolitonS50_K8_delta" + str(delt) + "_" + a + ".pdf", bbox_inches="tight", )
            plt.savefig("../plotDists/RobustSolitonS50_K8_delta" + str(delt) + "_" + a + ".svg")

            plt.close()

    S = 50
    for a in ["log", ""]:
        robust = IdealSolitonDistribution(S=S, seed=0)
        print(robust.pre_comp_dist)
        if a == "log":
            plt.semilogy([0.0] + robust.pre_comp_dist)
        else:
            plt.plot([0.0] + robust.pre_comp_dist)
        plt.ylabel("Probability")
        plt.xlabel("Degree")
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.show(block=False)
        plt.savefig("../plotDists/IdealSolitonS50_" + a + ".pdf", bbox_inches="tight")
        plt.savefig("../plotDists/IdealSolitonS50_" + a + ".svg")

        plt.close()


def onlineDist():
    for a in ["log", ""]:
        for eps in np.arange(0.01, 0.1, 0.01):
            dist = OnlineDistribution(eps=eps, seed=0)
            print(dist.pre_comp_dist)
            if a == "log":
                plt.semilogy([0.0] + dist.pre_comp_dist)
            else:
                plt.plot([0.0] + dist.pre_comp_dist)
            plt.ylabel("Probability")
            plt.xlabel("Degree")
            manager = plt.get_current_fig_manager()
            manager.resize(*manager.window.maxsize())
            plt.grid(True)
            plt.tight_layout()
            plt.show(block=False)
            plt.savefig("../plotDists/Online_eps" + str(eps) + "_" + a + ".pdf", bbox_inches="tight", )
            plt.savefig("../plotDists/Online_eps" + str(eps) + "_" + a + ".svg")
            plt.close()


def OnlineVergleich():
    for a in ["log", ""]:
        dist = OnlineDistribution(eps=0.01, seed=0)
        print("eps = 0.01")
        print(dist.pre_comp_dist[:30])
        if a == "log":
            plt.semilogy([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.01")
        else:
            plt.plot([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.01")

        dist = OnlineDistribution(eps=0.03, seed=0)
        print("eps = 0.03")
        print(dist.pre_comp_dist[:30])
        if a == "log":
            plt.semilogy([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.03")
        else:
            plt.plot([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.03")

        dist = OnlineDistribution(eps=0.06, seed=0)
        print("eps = 0.06")
        print(dist.pre_comp_dist[:30])
        if a == "log":
            plt.semilogy([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.06")
        else:
            plt.plot([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.06")

        plt.ylabel("Probability")
        plt.xlabel("Degree")
        plt.xticks(np.arange(0, 31, step=5))
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.show(block=False)
        plt.savefig("../plotDists/Vergleich_Online_30" + "_" + a + ".pdf", bbox_inches="tight")
        plt.savefig("../plotDists/Vergleich_Online_30" + "_" + a + ".svg")
        plt.close()


def raptorDist():
    for a in ["log", ""]:
        dist = RaptorDistribution(1)
        erg = [dist.deg(eps) for eps in [0, 10241, 491582, 712794, 831695, 948446, 1032189, 1048576]]
        print(erg)
        if a == "log":
            plt.semilogy([0.0] + erg)
        else:
            plt.plot([0.0] + erg)
        plt.ylabel("random number")
        plt.xlabel("Degree")
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.show(block=False)
        plt.savefig("../plotDists/raptor" + "_" + a + ".pdf", bbox_inches="tight")
        plt.savefig("../plotDists/raptor" + "_" + a + ".svg")
        plt.close()


def erlich_zielinski_robust_soliton_dist():
    for a in ["log", ""]:
        dist = ErlichZielinskiRobustSolitonDistribution(k=100, c=0.1, delta=0.05, seed=0)
        print(dist.pre_comp_dist)
        if a == "log":
            plt.semilogy([0.0] + dist.pre_comp_dist, label="log ErlichZlinski")
        else:
            plt.plot([0.0] + dist.pre_comp_dist, label="ErlichZlinski")
        """
        dist = OnlineDistribution(eps=0.03, seed=0)
        print("eps = 0.03")
        print(dist.pre_comp_dist[:30])
        if a == "log":
            plt.semilogy([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.03")
        else:
            plt.plot([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.03")

        dist = OnlineDistribution(eps=0.06, seed=0)
        print("eps = 0.06")
        print(dist.pre_comp_dist[:30])
        if a == "log":
            plt.semilogy([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.06")
        else:
            plt.plot([0.0] + dist.pre_comp_dist[:30], label="$\epsilon$ = 0.06")
    """
        plt.ylabel("Probability")
        plt.xlabel("Degree")
        plt.xticks(np.arange(0, 101, step=5))
        manager = plt.get_current_fig_manager()
        # manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.show(block=True)

        plt.close()


def raptorDistNeu():
    dist = [0.0, 10241 / 1048576, 481341 / 1048576, 221212 / 1048576, 118901 / 1048576, 0.0, 0.0, 0.0, 0.0, 0.0,
            116751 / 1048576, 83743 / 1048576, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16387 / 1048576, ]
    for a in ["log", ""]:
        print(dist)
        if a == "log":
            plt.semilogy(dist)
        else:
            plt.plot(dist)
        plt.ylabel("Wahrscheinlichkeit")
        plt.xlabel("Grad")
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.show(block=False)
        plt.savefig("../plotDists/raptor" + "_" + a + ".pdf", bbox_inches="tight")
        plt.savefig("../plotDists/raptor" + "_" + a + ".svg")
        plt.close()


if __name__ == "__main__":
    # main()
    erlich_zielinski_robust_soliton_dist()
    vergleich()
