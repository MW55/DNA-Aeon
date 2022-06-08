import itertools
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
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
    name6 = "ALL_DNA_SINGLE_0906/o_out_simple.csv"

    df6 = pd.read_csv(name6, delimiter=",", engine="python",
                      usecols=["filename", "overhead", "number_of_chunks", "dec_input", "invalid_drop", "seed",
                               "result", "timeNeeded", "Kodierung", ], index_col=False, )

    fig, axs = plt.subplots(1, 1)
    axs.set_ylabel("Time to Decode")
    print(df6)
    df6.boxplot(column="overhead", rot=17.5, by=["result", "Kodierung"])
    plt.title("Expected outcome per overhead and encoding")

    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=1))
    plt.xlabel("")
    # get rid of the automatic 'Boxplot grouped by group_by_column_name' title
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("simplednarules/box_result_overhead_kodierung_neu.pdf",
                bbox_inches="tight", )
    plt.savefig("simplednarules/box_result_overhead_kodierung_neu.svg")

    df6["invalid_dropPercent"] = df6["invalid_drop"] / (
            df6["dec_input"] + df6["invalid_drop"]
    )
    df6["encodedPackets"] = (df6["overhead"] + 1.0) * df6["number_of_chunks"]
    print((df6["dec_input"] + df6["invalid_drop"]))
    print(df6["encodedPackets"])
    df6["Fehlerwahrscheinlichkeit"] = df6["invalid_drop"] / df6["encodedPackets"]

    df6.boxplot(column="error probability", by=["Kodierung"])
    plt.title("Average of dropped packets")
    plt.xlabel("")
    plt.suptitle("")
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=1))
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("simplednarules/box_fehlerwkeit_kodierung_neu.pdf", bbox_inches="tight", )
    plt.savefig("simplednarules/box_fehlerwkeit_kodierung_neu.svg")
    plt.close()

    mark = itertools.cycle(["o", "v", "P", "*", "h", "D", "X", "d", "|"])
    linest = itertools.cycle(["--", "-"])
    color = itertools.cycle(["m", "m", "r", "r", "g", "g", "b", "b", "y", "y", "k", "k"])
    fig, axs = plt.subplots(1, 1)
    for name, group in df6.groupby(["Kodierung", "result"]):
        print(group)
        name, ergebnis = name
        if not ergebnis:
            name = "_nolegend_"
        else:
            name = name.replace("_", " ")
        group.plot.kde(
            x="error probability",
            y="overhead",
            label=name,
            ax=axs,
            color=next(color),
            linestyle=next(linest),
        )
    data = {
        "Fehlerwahrscheinlichkeit": [0, 1, 2, 3, 4],
        "overhead": [-100, -200, -20, -120, -890],
    }

    pd.DataFrame.from_dict(data).plot.kde(
        x="error probability",
        y="overhead",
        label="",
        ax=axs,
        color="white",
        linestyle="-",
    )

    pd.DataFrame.from_dict(data).plot.kde(
        x="error probability",
        y="overhead",
        label="success",
        ax=axs,
        color="black",
        linestyle="-",
    )
    pd.DataFrame.from_dict(data).plot.kde(
        x="error probability",
        y="overhead",
        label="no success",
        ax=axs,
        color="black",
        linestyle="--",
    )
    plt.title("")
    plt.xlabel("")
    plt.ylabel("density")
    plt.suptitle("")
    plt.xlim(-0.05, 0.7)
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("simplednarules/density.pdf", bbox_inches="tight", )
    plt.savefig("simplednarules/density.svg")
    plt.close()

    #### ohne kde

    mark = itertools.cycle(["o", "v", "P", "*", "h", "D", "X", "d", "|"])
    linest = itertools.cycle(["--", "-"])
    color = itertools.cycle(["m", "m", "r", "r", "g", "g", "b", "b", "y", "y", "k", "k"])
    fig, axs = plt.subplots(1, 1)
    for name, group in df6.groupby(["Kodierung", "result"]):
        print(group)
        name, ergebnis = name
        if not ergebnis:
            name = "_nolegend_"
        else:
            name = name.replace("_", " ")
        group.groupby("overhead").count().number_of_chunks.plot(
            label=name, ax=axs, color=next(color), linestyle=next(linest)
        )
    axs.plot([0], [0], color="b", linestyle="-", label="success")
    axs.plot([0], [0], color="b", linestyle="--", label="no success")
    plt.title("")
    plt.xlabel("")
    plt.ylabel("density")
    plt.suptitle("")
    plt.legend()
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("simplednarules/densitsdfsafy.pdf", bbox_inches="tight", )
    plt.savefig("simplednarules/densitsadfafy.svg")
    plt.close()


plot()
