import json
from pathlib import Path

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Bitstream Vera Sans"
matplotlib.rcParams["mathtext.it"] = "Bitstream Vera Sans:italic"
matplotlib.rcParams["mathtext.bf"] = "Bitstream Vera Sans:bold"


def lin_fun(x, a, b):
    return np.multiply(a, x) + b


def quad_fun(x, a, b):
    return np.multiply(a, np.power(x, 2)) + b


def pow_three_fun(x, a, b):
    return np.multiply(a, np.power(x, 3)) + b


def pow_four_fun(x, a, b):
    return np.multiply(a, np.power(x, 4)) + b


def runtime():
    plt.rcParams.update({"font.size": 15})
    input_filen = "/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/runtime/rt-original-1k-01.json"
    x_data = []
    y_data = []
    pers = []
    orders = []

    with Path(input_filen).open as data_file:
        data = json.load(data_file)

    for per in data:
        per_x_data = []
        per_y_data = []

        for c in data[per]:
            per_x_data.append(float(c))
            per_y_data.append(float(data[per][c]) / 1000.0)
        pers.append(float(per))
        x_data.append(per_x_data)
        y_data.append(per_y_data)

    xfinal = []
    yfinal = []
    cols = []
    colors = mcolors.BASE_COLORS

    for i in range(len(x_data)):
        col, val = colors.popitem()
        if col == "w" or col == "k":
            col, val = colors.popitem()
            if col == "w" or col == "k":
                col, val = colors.popitem()
                cols.append(col)
        label = "% 6.3f" % pers[i]
        orders.append(np.argsort(x_data[i]))
        xfinal.append(np.array(x_data[i])[orders[i]])
        yfinal.append(np.array(y_data[i])[orders[i]])
        plt.plot(xfinal[i], yfinal[i], "o", label="p=" + label, color=col)
        # start = 5
        # optimizedParameters, pcov = opt.curve_fit(LinFun, xfinal[i][start:], yfinal[i][start:])
        # plt.plot(xfinal[i][start:], LinFun(xfinal[i][start:], *optimizedParameters), '--', color=col, label='O(n)')
        print(xfinal)
        print(yfinal)

    # start2 = 7
    # optimizedParameters2, pcov = opt.curve_fit(LinFun, xfinal[1][start2:], yfinal[1][start2:])
    # plt.plot(xfinal[1][start2:], LinFun(xfinal[1][start2:], *optimizedParameters2),'--', color='m', label='O(n)')

    plt.legend()
    plt.ylabel("avg runtime (s) to decode 1k samples")
    plt.xlabel("code length n")
    plt.grid()
    plt.show()


def runtime_comparison():
    plt.rcParams.update({"font.size": 14})
    input_filen = ""
    input_filen2 = ""
    x_data = []
    y_data = []
    pers = []
    x_data2 = []
    y_data2 = []
    pers2 = 0.0
    orders = []

    with Path(input_filen).open() as data_file:
        data = json.load(data_file)
    with Path(input_filen2).open() as data_file2:
        data2 = json.load(data_file2)

    for per in data:
        per_x_data = []
        per_y_data = []

        for c in data[per]:
            per_x_data.append(float(c))
            per_y_data.append(float(data[per][c]) / 1000.0)
        pers.append(float(per))
        x_data.append(per_x_data)
        y_data.append(per_y_data)

    for per in data2:

        for c in data2[per]:
            x_data2.append(float(c))
            y_data2.append(float(data2[per][c]) / 1000.0)
        pers2 = float(per)

    xfinal = []
    yfinal = []
    cols = []
    colors = mcolors.BASE_COLORS

    for i in range(len(x_data)):
        col, val = colors.popitem()
        if col == "w" or col == "k":
            col, val = colors.popitem()
            if col == "w" or col == "k":
                col, val = colors.popitem()
                cols.append(col)
        label = "%2.2f" % pers[i]
        orders.append(np.argsort(x_data[i]))
        xfinal.append(np.array(x_data[i])[orders[i]])
        yfinal.append(np.array(y_data[i])[orders[i]])
        plt.plot(xfinal[i], yfinal[i], "o", label="UFH, p=" + label, color=col)
        start = 0
        optimized_parameters, pcov = opt.curve_fit(lin_fun, xfinal[i][start:], yfinal[i][start:])
        plt.plot(xfinal[i][start:], lin_fun(xfinal[i][start:], *optimized_parameters), "--", color=col, label=r"$O(n)$")
        print(xfinal)
        print(yfinal)

    # general qlpd decoder data
    label = "%2.2f" % pers2
    orders2 = np.argsort(x_data2)
    xfinal2 = np.array(x_data2)[orders2]
    yfinal2 = np.array(y_data2)[orders2]
    plt.plot(xfinal2, yfinal2, "d", label="GD, p=" + label, color="green")
    start = 0
    optimized_parameters2, pcov = opt.curve_fit(quad_fun, xfinal2[start:], yfinal2[start:])
    plt.plot(
        xfinal2[start:], quad_fun(xfinal2[start:], *optimized_parameters2), "--", color="green", label=r"$O(n^{2})$"
    )

    plt.legend()
    plt.ylabel("avg runtime (s) to decode $10^3$ samples")
    plt.xlabel("code length n")
    plt.grid()
    plt.show()


runtime_comparison()
