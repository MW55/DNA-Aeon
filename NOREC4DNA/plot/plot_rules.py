import os
import matplotlib
import pandas as pd
from cycler import cycler
import matplotlib.pyplot as plt

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


def custom_round(x, base=5):
    return int(base * round(float(x) / base))


"""
['seaborn-dark', 'classic', 'seaborn', 'fivethirtyeight',
'seaborn-dark-palette', 'dark_background', 'seaborn-poster',
'ggplot', 'seaborn-whitegrid', 'seaborn-colorblind', 'grayscale',
'seaborn-talk', 'seaborn-ticks', 'seaborn-paper', 'seaborn-deep',
'seaborn-darkgrid', 'seaborn-notebook', 'bmh', 'seaborn-bright',
'seaborn-muted', 'seaborn-white', 'seaborn-pastel']
"""

fname = "D:/Users/thejanky/Desktop/Uni/WS 17-18/CSV/recent/ALL_DNA_SINGLE_0906/ALL_DNA_SIMPLE_025_schritte/dna_abc.csv"
df = pd.read_csv(fname, delimiter=",", engine="python", index_col=False,
                 usecols=["Algorithm", "A_Permutation", "T_Permutation", "C_Permutation", "G_Permutation",
                          "dinucleotid_Runs", "Homopolymers", "GC_Content", "Trinucleotid_Runs", "Did_Drop", ], )

if not os.path.isdir("pdfs"):
    os.makedirs("pdfs")
font = {"family": "normal", "size": 14}
axis_font = {"size": "14"}
matplotlib.rc("font", **font)

plt.rc("axes",
       prop_cycle=(cycler("color",
                          ["m", "m", "r", "r", "r", "c", "c", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y", "y",
                           "y", "k", "k", "k", "k", ], ) + cycler("linestyle",
                                                                  ["-", "--", "-", "--", ":", "-", "--", "-", "--", ":",
                                                                   "-.", "-", "--", ":", "-.", "-", "--", ":", "-.",
                                                                   "-", "--", ":", "-.", ], )), )

plt.rc("grid", c="0.5", ls=":", lw=1)
plt.grid(True)
plt.show(block=False)

for name, group in df.groupby(["Algorithm"]):
    group.boxplot(by="Did_Drop", showfliers=False)
    manager = plt.get_current_fig_manager()
    print(manager.window.maxsize())
    print(*manager.window.maxsize())
    print(type(manager.window.maxsize()))
    manager.resize(*(3840, 2160))
    plt.tight_layout()
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("abc" + name + ".pdf", bbox_inches="tight")
    plt.savefig("abc" + name + ".svg")
    plt.close()
input("...")
