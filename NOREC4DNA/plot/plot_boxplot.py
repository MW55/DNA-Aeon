import pandas as pd
from math import floor
import matplotlib.pyplot as plt

# Get current size
fig_size = plt.rcParams["figure.figsize"]
fig_dpi = plt.rcParams["figure.dpi"]
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

fname = "../../../../CSV/brauchbar/merge_postprocessed.csv"  # union.csv_comb'
df = pd.read_csv(fname, delimiter=",", engine="python", index_col=False)
fname2 = "../../../../CSV/brauchbar/merge_postprocessed.csv"  # union.csv_comb'
df_einzeln = pd.read_csv(fname2, delimiter=",", engine="python", index_col=False)
print(df.keys())
days = str(floor(df["timeNeeded"].sum() / 60 / 60 / 24))
hours = str(floor(df["timeNeeded"].sum() % (60 * 60 * 24) / 60 / 60))
minutes = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) / 60)
secs = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) % 60)
print(str(days) + ":" + str(hours) + ":" + str(minutes) + ":" + str(secs))
df.sort_values(by=["codecName", "numberOfEncodedPackets", "number_of_chunks"], inplace=True)  # , 'droprate'

tmp = df.plot(subplots=True, grid=True, kind="box", by="numberOfEncodedPackets", title=" ")
plt.grid(True)
plt.show(block=True)
