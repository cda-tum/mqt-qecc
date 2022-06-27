import json
import math

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


def LinFun(x, a, b):
    return np.multiply(a, x) + b


def QuadFun(x, a, b):
    return np.multiply(a, np.power(x, 2)) + b


def PowThreeFun(x, a, b):
    return np.multiply(a, np.power(x, 3)) + b


def runtime():
    inputFilen = '/home/luca/Documents/codeRepos/qunionfind/scripts/raw-final27-06-2022.json'
    fig, ax = plt.subplots()
    colors = mcolors.BASE_COLORS
    xData = []
    yData = []
    pers = []
    orders = []
    with open(inputFilen) as data_file:
        data = json.load(data_file)

    for per in data:
        perXData = []
        perYData = []

        for c in data[per]:
            if float(data[per][c]) != 0.0:
                perXData.append(float(c))
                perYData.append(float(data[per][c]))
        pers.append(float(per))
        xData.append(perXData)
        yData.append(perYData)

    for i in range(len(xData)):
        col, val = colors.popitem()
        if (col == 'w' or col == 'k'):
            col, val = colors.popitem()
            if (col == 'w' or col == 'k'):
                col, val = colors.popitem()
        label = '% 6.3f' % pers[i]
        orders.append(np.argsort(xData[i]))
        xData[i] = np.array(xData[i])[orders[i]]
        yData[i] = np.array(yData[i])[orders[i]]
        print(xData[i])
        print(yData[i])
        ax.plot(xData[i], yData[i], 'o', label='PER=' + label, color='b')
        print(xData[i])
        print(yData[i])
        linx = xData[i].copy()
        liny = yData[i].copy()
        optimizedParameters, pcov = opt.curve_fit(LinFun, linx, liny)
        optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, linx, liny)
        print(str(math.ceil(optimizedParameters[1])) + 'x+' + str(math.ceil(optimizedParameters[0])))
        ax.plot(linx, LinFun(linx, *optimizedParameters), color='b', label='O(n)')
        ax.plot(linx, QuadFun(linx, *optimizedParameters2), color='r', label='O(n)')
        # ax.plot(xData[i], np.power(xData[i],3), color='r', label='O(n^3)')

    # fits

    ax.set_xlabel('code size')
    ax.set_ylabel('Avg runtime(ms)')
    ax.legend()
    plt.show()


runtime()
