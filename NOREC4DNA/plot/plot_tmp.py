import os
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
font = {"family": "normal", "size": 14}
axis_font = {"size": "14"}
matplotlib.rc("font", **font)
plt.rcParams["svg.fonttype"] = "none"
# Prints: [8.0, 6.0]
print("Current size:", fig_size)
print("Current dpi:", fig_dpi)
# Set figure width to 12 and height to 9
fig_size[0] = 16
fig_size[1] = 9
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


def plot(fname):
    name4 = "../../../CSV/recent/lt_ideal_700_840_vorkommen_je_chunk.csv"
    name5 = "../../../CSV/recent/lt_robust_700_840_vorkommen_je_chunk.csv"
    name1 = "../../../CSV/recent/raptor_800_801_vorkommen_je_chunk.csv"
    name2 = "../../../CSV/recent/raptor_700_701_vorkommen_je_chunk.csv"
    name3 = "../../../CSV/recent/online_388_485_vorkommen_je_chunk.csv"

    df1 = pd.read_csv(name1, delimiter=",", engine="python", index_col=False)
    df2 = pd.read_csv(name2, delimiter=",", engine="python", index_col=False)
    df3 = pd.read_csv(name3, delimiter=",", engine="python", index_col=False)
    df4 = pd.read_csv(name4, delimiter=",", engine="python", index_col=False)
    df5 = pd.read_csv(name5, delimiter=",", engine="python", index_col=False)
    df = df1.join(df2).join(df3).join(df4).join(df5)

    if not os.path.isdir("pdfs1"):
        os.makedirs("pdfs1")
    print(df)
    tmp = df.plot(subplots=True, grid=True, kind="box", sharey=True, widths=0.5)
    plt.grid(True)
    plt.show(block=False)
    manager = plt.get_current_fig_manager()
    print(manager.window.maxsize())
    print(*manager.window.maxsize())
    print(type(manager.window.maxsize()))
    manager.resize(*manager.window.maxsize())
    plt.tight_layout()
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("pdfs1/box_amount_per_chunk_" + fname.split(".")[0] + ".pdf", bbox_inches="tight", )
    plt.savefig("pdfs1/box_amount_per_chunk_" + fname.split(".")[0] + ".svg")
    plt.close()


name1 = "raptor_800_801_vorkommen_je_chunk.csv"  # union.csv_comb'
name2 = "lt_robust_700_935_vorkommen_je_chunk.csv"
name3 = "lt_ideal_700_935_vorkommen_je_chunk.csv"
name4 = "raptor_700_701_vorkommen_je_chunk.csv"
name5 = "online_388_485_vorkommen_je_chunk.csv"
plot("ALL")
