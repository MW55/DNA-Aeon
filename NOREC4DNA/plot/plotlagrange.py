import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 16
fig_size[1] = 9

plt.rcParams["svg.fonttype"] = "none"
fig_dpi = 500
plt.rcParams["figure.figsize"] = fig_size
plt.rcParams["figure.dpi"] = fig_dpi
font = {"family": "normal", "size": 24}
axis_font = {"size": "24"}
matplotlib.rc("font", **font)


def f(x):
    return (-(563 * pow(x, 9)) / 362880
            + (433 * pow(x, 8)) / 5760
            - (13421 * pow(x, 7)) / 8640
            + (51473 * pow(x, 6)) / 2880
            - (433411 * pow(x, 5)) / 3456
            + (3182497 * pow(x, 4)) / 5760
            - (8594527 * pow(x, 3)) / 5670
            + (1187297 * pow(x, 2)) / 480
            - (775529 * x) / 360 + 754)


xx = np.linspace(1.5, 9.5, 9)
x = np.linspace(1.0, 10.0, 10)
print(x)
print(xx)
jo = np.append(sorted(list(set(xx.tolist()) | set(x.tolist()))), [])
print(jo)
x2 = np.linspace(1.0, 10.0, 200)

y = f(jo)
p = lagrange(jo, y)
plt.plot(x, f(x), "o", label="support points", markersize=12)
plt.plot(xx, f(xx), "o", label="additional points", markersize=12)
plt.plot(x2, np.polyval(p, x2), label="Polynomial")
plt.xticks(np.arange(1.0, 11.0, 1))
plt.grid(True)
plt.tight_layout()
plt.legend(loc="upper left")


def g():
    return [f(a) for a in np.arange(1.0, 5.0, 0.0001)]


manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.show(block=False)
plt.savefig("Lagrange_eng.svg")
plt.savefig("Lagrange_eng.pdf", bbox_inches="tight")
plt.close()
