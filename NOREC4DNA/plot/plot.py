import os
import matplotlib
import pandas as pd
from math import floor
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

fname = "../../../CSV/csvs/merge.csv"  # union.csv_comb'
df = pd.read_csv(fname, delimiter=", ", engine="python", index_col=False)
fname2 = "../../../CSV/csvs/merge_einzeln.csv"  # union.csv_comb'
df_einzeln = pd.read_csv(fname2, delimiter=", ", engine="python", index_col=False)

print(str(df["timeNeeded"].sum()))
days = str(floor(df["timeNeeded"].sum() / 60 / 60 / 24))
hours = str(floor(df["timeNeeded"].sum() % (60 * 60 * 24) / 60 / 60))
minutes = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) / 60)
secs = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) % 60)

df.sort_values(
    by=["codecName", "numberOfEncodedPackets", "number_of_chunks"], inplace=True
)  # , 'droprate'

df["result"] = df["result"].astype(int)  # ([' result'])) #.cast('int')
df["codecName"] = df["codecName"].astype(str)
df["Overhead"] = (
        (df["numberOfEncodedPackets"] - df["number_of_chunks"]) - df["droppedCount"]
).apply(lambda x: custom_round(x, base=20))
df_einzeln.sort_values(
    by=["codecName", "numberOfEncodedPackets", "number_of_chunks"], inplace=True
)  # , 'droprate'
df_einzeln["result"] = df_einzeln["result"].astype(int)  # ([' result'])) #.cast('int')
df_einzeln["codecName"] = df_einzeln["codecName"].astype(str)
df_einzeln["Overhead"] = (
        (df_einzeln["numberOfEncodedPackets"] - df_einzeln["number_of_chunks"])
        - df_einzeln["droppedCount"]
).apply(lambda x: custom_round(x, base=20))

if not os.path.isdir("pdfs"):
    os.makedirs("pdfs")
font = {"family": "normal", "size": 14}
axis_font = {"size": "14"}
matplotlib.rc("font", **font)

# Nach codecName und number_of_chunks gruppiert
plt.rc(
    "axes",
    prop_cycle=(cycler("color",
                       ["m", "m", "r", "r", "r", "c", "c", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y", "y", "y",
                        "k", "k", "k", "k", ], )
                + cycler("linestyle", ["-", "--", "-", "--", ":", "-", "--", "-", "--", ":", "-.", "-", "--", ":", "-.",
                                       "-", "--", ":", "-.", "-", "--", ":", "-.", ], )),
)

plt.rc("grid", c="0.5", ls=":", lw=1)
plt.grid(True)
plt.show(block=False)

tmp1 = (
    df_einzeln.groupby(["droprate", "codecName"]).mean().unstack().plot(y="result", title="Erfolgsrate bzgl. droprate"))

manager = plt.get_current_fig_manager()
print(manager.window.maxsize())
print(*manager.window.maxsize())
print(type(manager.window.maxsize()))
manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
plt.tight_layout()
plt.title("Erfolg je Codec")
plt.ylabel("Erfolgsrate", **axis_font)
plt.grid(True)
plt.show(block=False)
plt.savefig("../pdfs/all_numchunks_" + ".pdf", bbox_inches="tight")
plt.savefig("../pdfs/all_numchunks_" + ".svg")
plt.close()

# Nur nach codecName Gruppiert
plt.rc("axes",
       prop_cycle=(cycler("color", ["r", "r", "r", "r", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y", "y", "y",
                                    "k", "k", "k", "k", ], ) + cycler("linestyle",
                                                                      ["-", "--", ":", "-.", "-", "--", ":", "-.",
                                                                       "-", "--", ":", "-.", "-", "--", ":", "-.",
                                                                       "-", "--", ":", "-.", ], )), )

tmp1 = (df.groupby(["droprate", "codecName"]).mean().unstack().plot(y="result", title="Erfolgsrate bzgl. droprate"))

manager = plt.get_current_fig_manager()
manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
plt.tight_layout()
plt.title("Erfolg je Codec")
plt.ylabel("Erfolgsrate", **axis_font)
plt.grid(True)
plt.show(block=False)
plt.savefig("../pdfs/all_" + ".pdf", bbox_inches="tight")
plt.savefig("../pdfs/all_" + ".svg")
plt.close()

plt.rc("axes",
       prop_cycle=(cycler("color", ["b", "k", "r", "g", "y"]) + cycler("linestyle", ["-", "-", "-", "-", "-"])), )

for name, group in df.groupby(["codecName", "number_of_chunks"]):
    name, noChu = name
    df1 = group.groupby(["droprate"]).mean()
    df2 = group.groupby(["Overhead"]).mean()

    tmp = df.loc[df["codecName"] == name].loc[df["number_of_chunks"] == noChu]
    m = tmp["number_of_chunks"].mean()
    name = str(name) + "_number_of_chunks=" + str(noChu)

    col = ["g" if res else "r" for res in tmp["result"]]
    tmp1 = tmp.plot(
        kind="scatter",
        x="numberOfEncodedPackets",
        y="timeNeeded",
        color=col,
        title=name.replace("_", " ") + " ( " + str(len(tmp)) + " Elemente )",
    )
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    plt.grid(True)
    tmp1.grid(True)
    plt.show(block=False)
    # plt.tight_layout()
    plt.savefig("../pdfs/scatter_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/scatter_" + name + ".svg")
    plt.close()

    col = ["g" if res else "r" for res in tmp["result"]]
    tmp1 = tmp.plot(kind="scatter", x="numberOfEncodedPackets", y="droprate", color=col,
                    title=name.replace("_", " ") + " ( " + str(len(tmp)) + " Elemente )", )
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    plt.grid(True)
    tmp1.grid(True)
    plt.show(block=False)
    # plt.tight_layout()
    plt.savefig("../pdfs/scatter_droprate_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/scatter_droprate_" + name + ".svg")
    plt.close()

    col = ["g" if res else "r" for res in tmp["result"]]
    tmp1 = tmp.plot(kind="scatter", x="numberOfEncodedPackets", y="droppedCount", color=col,
                    title=name.replace("_", " ") + " ( " + str(len(tmp)) + " Elemente )", )
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    plt.grid(True)
    tmp1.grid(True)
    plt.show(block=False)
    plt.savefig("../pdfs/scatter_droppedCount_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/scatter_droppedCount_" + name + ".svg")
    plt.close()

    tmp1 = df2.reset_index().plot(x="Overhead", y="result",
                                  title=name.replace("_", " ") + " ( " + str(len(tmp)) + " Elemente )", )
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    plt.grid(True)
    tmp1.grid(True)
    plt.show(block=False)
    plt.savefig("../pdfs/overhead_vs_res_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/overhead_vs_res_" + name + ".svg")
    plt.close()

    tmp1 = df1.reset_index().plot(x="droprate", y="result",
                                  title=name.replace("_", " ") + " ( " + str(len(tmp)) + " Elemente )", legend=False, )
    tmp1.set_ylabel("result")
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    tmp1.grid(True)
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("../pdfs/drop_vs_res_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/drop_vs_res_" + name + ".svg")
    plt.close()

    m = df1["number_of_chunks"].mean()
    df1 = df1.drop("number_of_chunks", 1)
    df1 = df1.drop("seed", 1)
    tmp = df1.plot(subplots=True, grid=True, use_index="number_of_chunks", title=name.replace("_", " "), legend=True, )
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("../pdfs/multi_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/multi_" + name + ".svg")
    plt.close()

    tmp = df1.plot(subplots=True, grid=True, kind="box", by="numberOfEncodedPackets",
                   title=name.replace("_", " "), )  # + " ( " + str(len(tmp)) + " Elemente )")
    manager = plt.get_current_fig_manager()
    manager.resize(*(3840, 2160))  # manager.resize(*manager.window.maxsize())
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("../pdfs/box_" + name + ".pdf", bbox_inches="tight")
    plt.savefig("../pdfs/box_" + name + ".svg")
    plt.close()

input("Press Enter to exit...")
