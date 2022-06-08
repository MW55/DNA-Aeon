import os
import matplotlib
import pandas as pd
from math import floor
from cycler import cycler
import matplotlib.pyplot as plt


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

fname = (
    "../../../../CSV/ALL_DNA_DOUBLE_RU10FAST/ALL_DNA_DOUBLE_RU10FAST/merge.csv"
)  # union.csv_comb'
df = pd.read_csv(
    fname, delimiter=",", engine="python", index_col=False, encoding="utf-16"
)

days = str(floor(df["timeNeeded"].sum() / 60 / 60 / 24))
hours = str(floor(df["timeNeeded"].sum() % (60 * 60 * 24) / 60 / 60))
minutes = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) / 60)
secs = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) % 60)
df.sort_values(by=["codecName", "overhead", "number_of_chunks", "dec_input", "invalid_drop", "timeNeeded", ],
               inplace=True, )  # , 'droprate'
print(df["number_of_chunks"])
df["result"] = df["result"].astype(int)
df["codecName"] = df["codecName"].astype(str)
df["number_of_chunks"] = df["number_of_chunks"].astype(int)
df["overhead"] = df["overhead"].astype(float).apply(lambda x: 1.0 + x)
df["numberOfEncodedPackets"] = (df["overhead"] * df["number_of_chunks"]).apply(lambda x: custom_round(x, base=1))

if not os.path.isdir("pdfs"):
    os.makedirs("pdfs")
font = {"family": "normal", "size": 14}
axis_font = {"size": "14"}
matplotlib.rc("font", **font)

##### Nach codecName und number_of_chunks Gruppiert #####
plt.rc("axes",
       prop_cycle=(cycler("color", ["m", "m", "r", "r", "r", "r", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y",
                                    "y", "y", "k", "k", "k", "k", ], ) + cycler("linestyle",
                                                                                ["-", "--", "-", "--", ":", "-.", "-",
                                                                                 "--", ":", "-.", "-", "--", ":", "-.",
                                                                                 "-", "--", ":", "-.", "-", "--", ":",
                                                                                 "-.", ], )), )

plt.rc("grid", c="0.5", ls=":", lw=1)

##### Nur nach codecName Gruppiert #####
tmp1 = (df.groupby(["invalid_drop", "codecName"]).mean().unstack().plot(y="result", title="Success rate per drop rate"))
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.show(block=False)
plt.close()

plt.rc("axes", prop_cycle=(cycler("color",
                                  ["m", "m", "r", "r", "r", "r", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y", "y",
                                   "y", "k", "k", "k", "k", ], ) + cycler("linestyle",
                                                                          ["-", "--", "-", "--", ":", "-.", "-", "--",
                                                                           ":", "-.", "-", "--", ":", "-.", "-", "--",
                                                                           ":", "-.", "-", "--", ":", "-.", ], )), )
fig, axs = plt.subplots(1, 1)
axs.set_xlim(-0.1, 1)
df["overhead"] = df["overhead"].astype(float).apply(lambda x: x - 1.0, 0)
axs.set_ylabel("Erfolgsrate")
for name, group in df.groupby(["codecName"]):
    df2 = group.groupby(["overhead"]).mean()
    tmp = df.loc[df["codecName"] == name]
    m = tmp["number_of_chunks"].mean()
    tmp1 = df2.reset_index()
    print(tmp1["result"])
    a = tmp1.plot(x="overhead", y="result", label=name.replace("_", " ").replace("eps", "\epsilon"), ax=axs,
                  legend=True, )
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.grid(True)
plt.tight_layout()
plt.show(block=False)
plt.savefig("../pdfs/overhead_percent_vs_res_" + name + ".pdf", bbox_inches="tight")
plt.savefig("../pdfs/overhead_percent_vs_res_" + name + ".svg")
plt.close()

input("Press Enter to exit...")
