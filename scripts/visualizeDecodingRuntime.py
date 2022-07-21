import json

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
    plt.rcParams.update({'font.size': 15})
    inputFilen = '/home/luca/Documents/codeRepos/qecc/scripts/in.json'
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
            avgRuntime = float(data[per][c])
            if avgRuntime != 0.0:  # skip 0 runtime
                perXData.append(float(c))
                perYData.append(avgRuntime)
        pers.append(float(per))
        xData.append(perXData)
        yData.append(perYData)

    xfinal =[]
    yfinal =[]
    colors = mcolors.BASE_COLORS

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
        plt.plot(xData[i], yData[i], 'o', label='PER=' + label, color='b')
        optimizedParameters, pcov = opt.curve_fit(LinFun, xData[i], yData[i])
        optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, xData[i], yData[i])
        plt.plot(xData[i], LinFun(xData[i], *optimizedParameters), color='b', label='O(n)')
        plt.plot(xData[i], QuadFun(xData[i], *optimizedParameters2), color='r', label='O(n^2)')
        print(str(optimizedParameters[0]) + 'x+' + str(optimizedParameters[1]))
        xfinal.append(np.array(xData[i])[np.argsort(xData[i])])
        yfinal.append(np.array(yData[i])[np.argsort(xData[i])])


    plt.legend()
    plt.grid()
    plt.show()
runtime()
