import itertools
import pandas as pd
import matplotlib.pyplot as plt

# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
plt.rcParams["svg.fonttype"] = "none"
print("Current size:", fig_size)
print("Current dpi:", fig_dpi)


def custom_round(x, base=5):
    return int(base * round(float(x) / base))


# plt.style.use('fivethirtyeight')
"""
['seaborn-dark', 'classic', 'seaborn', 'fivethirtyeight',
'seaborn-dark-palette', 'dark_background', 'seaborn-poster',
'ggplot', 'seaborn-whitegrid', 'seaborn-colorblind', 'grayscale',
'seaborn-talk', 'seaborn-ticks', 'seaborn-paper', 'seaborn-deep',
'seaborn-darkgrid', 'seaborn-notebook', 'bmh', 'seaborn-bright',
'seaborn-muted', 'seaborn-white', 'seaborn-pastel']
"""


def plot():
    name6 = "../../../../CSV/brauchbar/ALL_DNA_DOUBLE1006/o_out_simple.csv"

    df6 = pd.read_csv(name6, delimiter=",", engine="python",
                      usecols=["filename", "overhead", "number_of_chunks", "Decoder Input", "invalid_drop", "Seed",
                               "Ergebnis", "Sekunden", "Kodierung", "createdPackets", "DropRate", ], index_col=False, )

    fig, axs = plt.subplots(1, 1)
    axs.set_ylabel("Dauer Decode")
    print(df6)
    df6.boxplot(column="overhead", rot=17.5, by=["Ergebnis", "Kodierung"])
    plt.title("Erwarteter Ausgang je Overhead und Kodierung")
    plt.xlabel("")
    # get rid of the automatic 'Boxplot grouped by group_by_column_name' title
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../TeX/Bilder/doublednarules/box_result_overhead_kodierung_neu.pdf", bbox_inches="tight", )
    plt.savefig("../../../TeX/Bilder/doublednarules/box_result_overhead_kodierung_neu.svg")

    df6["invalid_dropPercent"] = df6["invalid_drop"] / (df6["Decoder Input"] + df6["invalid_drop"])
    df6["encodedPackets"] = (df6["overhead"] + 1.0) * df6["number_of_chunks"]
    print((df6["Decoder Input"] + df6["invalid_drop"]))
    print(df6["encodedPackets"])
    df6["Fehlerwahrscheinlichkeit"] = df6["invalid_drop"] / df6["createdPackets"]

    df6.boxplot(column="Fehlerwahrscheinlichkeit", rot=17.5, by=["Kodierung"])
    plt.title("Vom Simulator durchschnittlich verworfene Pakete AAA")
    plt.xlabel("")
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../TeX/Bilder/doublednarules/box_fehlerwkeit_kodierung_neu.pdf", bbox_inches="tight", )
    plt.savefig("../../../TeX/Bilder/doublednarules/box_fehlerwkeit_kodierung_neu.svg")
    plt.close()

    mark = itertools.cycle(["o", "v", "P", "*", "h", "D", "X", "d", "|"])
    linest = itertools.cycle(["--", "-"])
    color = itertools.cycle(["m", "m", "r", "r", "g", "g", "b", "b", "y", "y", "k", "k"])
    fig, axs = plt.subplots(1, 1)
    for name, group in df6.groupby(["Kodierung", "Ergebnis"]):
        print(group)
        name, ergebnis = name
        if not ergebnis:
            name = "_nolegend_"
        else:
            name = name.replace("_", " ")
        group.plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label=name, ax=axs, color=next(color),
                       linestyle=next(linest), )
    data = {"Fehlerwahrscheinlichkeit": [0, 1, 2, 3, 4], "overhead": [-100, -200, -20, -120, -890], }

    pd.DataFrame.from_dict(data).plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label="", ax=axs, color="white",
                                          linestyle="-", )

    pd.DataFrame.from_dict(data).plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label="Erfolgreich", ax=axs,
                                          color="black", linestyle="-", )
    pd.DataFrame.from_dict(data).plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label="nicht Erfolgreich", ax=axs,
                                          color="black", linestyle="--", )
    plt.title("")
    plt.xlabel("")
    plt.ylabel("Dichte")
    plt.suptitle("")
    plt.xlim(-0.05, 0.6)
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../TeX/Bilder/doublednarules/density.pdf", bbox_inches="tight")
    plt.savefig("../../../TeX/Bilder/doublednarules/density.svg")
    plt.close()


plot()
