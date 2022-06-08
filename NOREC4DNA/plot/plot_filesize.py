import itertools
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["savefig.dpi"] = 10000
# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
fig_size[0] = 16
fig_size[1] = 9

plt.rcParams["svg.fonttype"] = "none"
fig_dpi = 500
plt.rcParams["figure.figsize"] = fig_size
plt.rcParams["figure.dpi"] = fig_dpi
plt.rcParams["svg.fonttype"] = "none"
# Prints: [8.0, 6.0]
print("Current size:", fig_size)
print("Current dpi:", fig_dpi)
msize = 10

font = {"family": "normal", "size": 38}
axis_font = {"size": "38"}
matplotlib.rc("font", **font)
"""
['seaborn-dark', 'classic', 'seaborn', 'fivethirtyeight',
'seaborn-dark-palette', 'dark_background', 'seaborn-poster',
'ggplot', 'seaborn-whitegrid', 'seaborn-colorblind', 'grayscale',
'seaborn-talk', 'seaborn-ticks', 'seaborn-paper', 'seaborn-deep',
'seaborn-darkgrid', 'seaborn-notebook', 'bmh', 'seaborn-bright',
'seaborn-muted', 'seaborn-white', 'seaborn-pastel']
"""


def plot_encode():
    name6 = "../../../../CSV/recent/DATEIGROESSE_Zeit/orange_out_new_new.csv"
    axis_font = {"size": "38"}

    df6 = pd.read_csv(name6, delimiter=",", engine="python",
                      usecols=["Art", "number_of_chunks", "Sekunden", "Datei", "Modus"], index_col=False, )
    df6["Sekunden"] = 1.0 * df6["Sekunden"] / 10.0
    fig, axs = plt.subplots(1, 1)

    for name, group in df6.groupby(["Art"]):
        fig, axs = plt.subplots(1, 1)
        axs.set_ylabel("Sekunden", **axis_font)
        mark = itertools.cycle(["o", "v", "^", ">", "1", "2", "3", "4", "8", "P", "*", "h", "D", "X", "d", "|", ])
        linest = itertools.cycle(["-", ":"])
        color = itertools.cycle(["m", "m", "r", "r", "c", "c", "g", "g", "b", "b", "y", "y", "k", "k"])
        for name1, group1 in group.groupby(["Datei", "Modus"]):
            datei, modus = name1
            tmp = group1.groupby(["number_of_chunks"]).mean()
            yyy = tmp.reset_index()
            if modus == "ohne Speichern":
                datei = "_nolegend_"
            trenn = ""
            modus = ""
            yyy.plot(x="number_of_chunks", y="Sekunden", label=datei + trenn + modus, ax=axs, color=next(color),
                     marker=next(mark), markersize=msize, linestyle=next(linest), linewidth=3.0, )
        axs.plot([], [], label="   ", color="white", linestyle="-")
        axs.plot([], [], color="black", linestyle="-", label="ohne IO")
        axs.plot([], [], color="black", linestyle="--", label="ohne Speichern")
        plt.title("Encode: " + name)
        plt.xlabel("Anzahl Chunks", **axis_font)
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.legend(loc=2)
        plt.tight_layout()
        plt.show(block=False)
        fig = plt.gcf()
        plt.savefig("size_Zeit/Encode/chunks_vs_time_" + name + "_bigger.pdf", bbox_inches="tight", )
        plt.savefig("size_Zeit/Encode/chunks_vs_time_" + name + "_bigger.svg")
    df6 = pd.read_csv(name6, delimiter=",", engine="python", usecols=["Art", "Kilobytes", "Sekunden", "Datei", "Modus"],
                      index_col=False, )
    df6["Sekunden"] = 1.0 * df6["Sekunden"] / 10.0
    plt.close()
    fig, axs = plt.subplots(1, 1)
    linest = itertools.cycle(["-", ":"])
    color = itertools.cycle(["m", "m", "r", "r", "c", "c", "g", "g", "b", "b", "y", "y", "k", "k"])
    for name, group in df6.groupby(["Art"]):

        for name1, group1 in group.groupby(["Modus"]):
            modus = name1
            tmp = group1.groupby(["Kilobytes"]).mean()
            y = tmp.reset_index()
            trenn = " - "
            y.plot(x="Kilobytes", y="Seconds", label=name + trenn + modus, ax=axs, legend=True,
                   linestyle=next(linest), linewidth=3.0, color=next(color), )
            plt.title("Time needed compared to the filesize")
    plt.ylabel("Seconds", **axis_font)
    plt.xlabel("Kilobytes", **axis_font)
    plt.legend(loc=2)
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../../Masterarbeit/TeX/Bilder/DATEIGROESSE_Zeit/Encode/Kilobytes_vs_time_" + "_bigger.pdf",
                bbox_inches="tight", )
    plt.savefig("../../../../Masterarbeit/TeX/Bilder/DATEIGROESSE_Zeit/Encode/Kilobytes_vs_time_" + "_bigger.svg")
    plt.close()

    df6 = pd.read_csv(name6, delimiter=",", engine="python",
                      usecols=["Art", "Kilobytes", "Sekunden", "Datei", "Modus", "number_of_chunks"], index_col=False, )
    df6["Sekunden"] = 1.0 * df6["Sekunden"] / 10.0
    fig, axs = plt.subplots(1, 1)
    linest = itertools.cycle(["-", ":"])

    mark = itertools.cycle(["o", "v", "^", ">", "1", "2", "3", "4", "8", "P", "*", "h", "D", "X", "d", "|"])
    color = itertools.cycle(["m", "m", "r", "r", "c", "c", "g", "g", "b", "b", "y", "y", "k", "k"])
    for name, group in df6.groupby(["Art"]):
        fig, axs = plt.subplots(1, 1)
        for name1, group1 in group.groupby(["number_of_chunks"]):
            tmp = group1.groupby(["Kilobytes"]).mean()
            y = tmp.reset_index()
            y.plot(x="Kilobytes", y="Seconds", label=name + " with " + str(int(name1)) + " Chunks", ax=axs, legend=True,
                   linestyle=next(linest), linewidth=3.0, color=next(color), marker=next(mark), markersize=msize, )

        plt.title("Time needed compared to filesize")
        plt.ylabel("Seconds", **axis_font)
        plt.xlabel("Kilobytes", **axis_font)
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.show(block=False)
        plt.savefig("size_Zeit//Encode/Kilobytes_vs_time_packets_" + name + "_bigger.pdf", bbox_inches="tight", )
        plt.savefig("size_Zeit/Encode/Kilobytes_vs_time_packets_" + name + "_bigger.svg")
        plt.close()


def plot_decode():
    name6 = "size_Zeit/orange_out_decode_ohneapprox_new.csv"
    axis_font = {"size": "38"}

    df6 = pd.read_csv(name6, delimiter=",", engine="python", index_col=False)
    df6["Sekunden"] = 1.0 * df6["Sekunden"] / 10.0
    fig, axs = plt.subplots(1, 1)
    for name, group in df6.groupby(["Art"]):
        fig, axs = plt.subplots(1, 1)
        axs.set_ylabel("Sekunden", **axis_font)
        mark = itertools.cycle(["o", "v", "P", "*", "h", "D", "X", "d", "|"])
        linest = itertools.cycle(["-", "-"])
        color = itertools.cycle(["m", "r", "c", "g", "b", "y", "k"])
        for name1, group1 in group.groupby(["Datei", "Modus"]):
            datei, modus = name1
            tmp = group1.groupby(["number_of_chunks"]).mean()
            yyy = tmp.reset_index()
            if modus == "ohne Speichern":
                datei = "_nolegend_"
            trenn = ""
            modus = ""
            yyy.plot(x="number_of_chunks", y="seconds", label=datei + trenn + modus, ax=axs, color=next(color),
                     marker=next(mark), markersize=msize, linestyle=next(linest), linewidth=3.0, )
        plt.title("Decoder: " + name)
        plt.xlabel("#Chunks", **axis_font)
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.legend(fontsize=16)
        plt.tight_layout()
        plt.show(block=False)
        fig = plt.gcf()
        plt.savefig("size_time/Decode/chunks_vs_time_" + name + "_bigger.pdf", bbox_inches="tight", )
        plt.savefig("size_time/Decode/chunks_vs_time_" + name + "_bigger.svg")
        plt.close()

        fig, axs = plt.subplots(1, 1)
        axs.set_ylabel("seconds", **axis_font)
        mark = itertools.cycle(["o", "v", "P", "*", "h", "D", "X", "d", "|"])
        linest = itertools.cycle(["-", "-"])
        color = itertools.cycle(["m", "r", "c", "g", "b", "y", "k"])
        for name1, group1 in group.groupby(["file", "mode"]):
            datei, modus = name1
            tmp = group1.groupby(["number_of_packets"]).mean()
            yyy = tmp.reset_index()
            trenn = " - "
            if modus == "ohne Speichern" or modus == "ohne IO":
                trenn = ""
                modus = ""
            yyy.plot(x="number_of_packets", y="seconds", label=datei + trenn + modus, ax=axs, color=next(color),
                     marker=next(mark), markersize=msize, linestyle=next(linest), linewidth=3.0, )
        plt.title("Decoder: " + name)
        plt.xlabel("#Pakets", **axis_font)
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        font = {"family": "normal", "size": 38}
        axis_font = {"size": "38"}
        matplotlib.rc("font", **font)
        plt.legend(fontsize=25)
        plt.tight_layout()
        plt.show(block=False)
        fig = plt.gcf()
        plt.savefig("size_time/Decode/big_packets_vs_time_" + name + "_bigger.pdf", bbox_inches="tight", )
        plt.savefig("size_time/Decode/big_packets_vs_time_" + name + "_bigger.svg")
        plt.close()


# plot_encode()
plot_decode()
