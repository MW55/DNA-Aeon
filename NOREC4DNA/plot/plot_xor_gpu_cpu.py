import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
plt.rcParams["svg.fonttype"] = "none"
# Prints: [8.0, 6.0]
print("Current size:", fig_size)
print("Current dpi:", fig_dpi)


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


def plot():
    name6 = "../cpu_gpu.csv"
    font = {"family": "normal", "size": 14}
    axis_font = {"size": "14"}
    matplotlib.rc("font", **font)

    df6 = pd.read_csv(name6, delimiter=",", engine="python",
                      usecols=["Name", "Anzahl Chunks", "Sekunden", "Datei", "Art"], index_col=False, )

    fig, axs = plt.subplots(1, 1)
    for name, group in df6.groupby(["Art"]):
        print(group.groupby(["Anzahl Chunks"]).mean().reset_index())
        group.groupby(["Anzahl Chunks"]).mean().reset_index().plot(x="#Chunks", y="Seconds", label=name, ax=axs)
    plt.title("")
    plt.xlabel("#Chunks", **axis_font)
    plt.ylabel("Seconds", **axis_font)
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("gpu_cpu_xor/in_framework_gpu_cpu.pdf", bbox_inches="tight", )
    plt.savefig("gpu_cpu_xor/in_framework_gpu_cpu.svg")
    plt.close()


plot()
