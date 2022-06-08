import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

# Get current size
# plt.style.use('seaborn-whitegrid')
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


def plot(column=True):
    name6 = "../../../../CSV/brauchbar/GPU_CPU_XOR/good.csv"
    font = {"family": "normal", "size": 18}
    axis_font = {"size": "18"}
    matplotlib.rc("font", **font)
    x_achse = "Count " + ("columns" if column else "rows")
    df6 = pd.read_csv(name6, delimiter=",", engine="python", usecols=["Method", x_achse, "Seconds"],
                      index_col=False, )  # "Anzahl Spalten",

    plt.title("Comparison of the XOR reductions")
    plt.xlabel("")
    plt.suptitle("")

    for name, group in df6.groupby(["Methode"]):
        group = group.groupby([x_achse]).mean()
        if name != "CUDA_SIMPLE":
            # hamming
            # blackman
            # bartlett
            # parzen
            # barthannc
            plt.plot(group, label=name)
            print(group.rolling(window=15, win_type="blackman", closed="neither", min_periods=1).mean(std=1))

    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.ylabel("Seconds", **axis_font)
    plt.xlabel(x_achse, **axis_font)
    plt.legend()
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("../../../../CSV/brauchbar/GPU_CPU_XOR/" + str(x_achse.replace(" ", "_")) + "_vs_zeit.pdf",
                bbox_inches="tight", )
    plt.savefig("../../../../CSV/brauchbar/GPU_CPU_XOR/" + str(x_achse.replace(" ", "_")) + "_vs_zeit.svg")
    plt.close()


plot()
