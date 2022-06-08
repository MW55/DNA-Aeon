import itertools
import pandas as pd
import matplotlib.pyplot as plt

# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
# Prints: [8.0, 6.0]
plt.rcParams["svg.fonttype"] = "none"
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
    name6 = "../../../../CSV/brauchbar/ALL_DNA_VERGLEICH/Orange_out_new.csv"

    df6 = pd.read_csv(name6, delimiter=",", engine="python",
                      usecols=["filename", "overhead", "number_of_chunks", "Decoder Input", "invalid_drop",
                               "Seed", "Ergebnis", "Sekunden", "createdPackets", "DropRate", "Kodierung", "Versuch",
                               "group", ], index_col=False, )

    fig, axs = plt.subplots(1, 1)
    axs.set_ylabel("Dauer Decode")
    print(df6)
    df6.boxplot(column="overhead", rot=17.5, by=["Ergebnis", "group"])
    plt.title("Erwarteter Ausgang je Overhead und Kodierung")
    plt.xlabel("")
    # get rid of the automatic 'Boxplot grouped by group_by_column_name' title
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../../Masterarbeit/TeX/Bilder/alldnarules/box_result_overhead_kodierung_neu.pdf",
                bbox_inches="tight", )
    plt.savefig("../../../../Masterarbeit/TeX/Bilder/alldnarules/box_result_overhead_kodierung_neu.svg")

    df6["invalid_dropPercent"] = df6["invalid_drop"] / (
            df6["Decoder Input"] + df6["invalid_drop"]
    )
    df6["encodedPackets"] = (df6["overhead"] + 1.0) * df6["number_of_chunks"]
    print((df6["Decoder Input"] + df6["invalid_drop"]))
    print(df6["encodedPackets"])
    df6["Fehlerwahrscheinlichkeit"] = df6["invalid_drop"] / df6["createdPackets"]

    df6.boxplot(column="Fehlerwahrscheinlichkeit", rot=17.5, by=["group"])
    plt.title("Vom Simulator durchschnittlich verworfene Pakete AAA")
    plt.xlabel("")
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../../Masterarbeit/TeX/Bilder/alldnarules/box_fehlerwkeit_kodierung_neu.pdf",
                bbox_inches="tight", )
    plt.savefig("../../../../Masterarbeit/TeX/Bilder/alldnarules/box_fehlerwkeit_kodierung_neu.svg")
    plt.close()

    for name, df6 in df6.groupby(["Kodierung"]):
        linest = itertools.cycle(["--", "-"])
        color = itertools.cycle(["m", "m", "r", "r", "g", "g", "b", "b", "y", "y", "k", "k", "c", "c"])
        fig, axs = plt.subplots(1, 1)
        for gname, group in df6.groupby(["group"]):
            for ergebnis, group in group.groupby(["Ergebnis"]):
                print(group)
                if not ergebnis:
                    gname1 = "_nolegend_"
                else:
                    gname1 = (gname.replace("_", " ")
                              .replace("LT Ideal", "")
                              .replace("LT Robust", "")
                              .replace("Online", "")
                              .replace("Raptor", "")
                              .replace("Double", "zweifach")
                              .replace("Scale", "skaliert")
                              .replace("Single", "einfach"))
                group.plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label=gname1, ax=axs, color=next(color),
                               linestyle=next(linest), )
        data = {"Fehlerwahrscheinlichkeit": [0, 1, 2, 3, 4],
                "overhead": [-100, -200, -20, -120, -890], }

        pd.DataFrame.from_dict(data).plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label="", ax=axs,
                                              color="white", linestyle="-", )

        pd.DataFrame.from_dict(data).plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label="Erfolgreich", ax=axs,
                                              color="black", linestyle="-", )
        pd.DataFrame.from_dict(data).plot.kde(x="Fehlerwahrscheinlichkeit", y="overhead", label="nicht Erfolgreich",
                                              ax=axs, color="black", linestyle="--", )
        plt.title(name)
        plt.xlabel("Overhead")
        plt.ylabel("Dichte")
        plt.suptitle("")
        plt.xlim(-0.05, 0.7)
        manager = plt.get_current_fig_manager()
        manager.resize(*manager.window.maxsize())
        plt.grid(True)
        plt.tight_layout()
        plt.show(block=False)
        plt.savefig("../../../../Masterarbeit/TeX/Bilder/alldnarules/density_" + name.replace(" ", "_") + ".pdf",
                    bbox_inches="tight", )
        plt.savefig("../../../../Masterarbeit/TeX/Bilder/alldnarules/density_" + name.replace(" ", "_") + ".svg")
        plt.close()


def plot2():
    # old = DNA_random_bytes_0_0.0_sim2018-07-12_20-36
    fname = ("../../../../Masterarbeit/Code/DNA_random_bytes_0_0.0_sim2018-07-13_01-10.csv")
    df6 = pd.read_csv(fname, delimiter=",", engine="python", index_col=False)

    fig, axs = plt.subplots(1, 1)
    df6.boxplot(column="Overall_Dropchance")  # , by=['Ergebnis', 'group'])
    plt.title("Verteilung der berechneten Dropchance je Paket")
    plt.xlabel("")
    # get rid of the automatic 'Boxplot grouped by group_by_column_name' title
    plt.suptitle("")
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig("../../../../Masterarbeit/Code/random_bytes.pdf", bbox_inches="tight")
    plt.savefig("../../../../Masterarbeit/Code/random_bytes.svg")


plot()
# plot2()
