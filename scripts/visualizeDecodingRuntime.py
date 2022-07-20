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
    inputFilen = '/home/luca/Documents/codeRepos/qecc/scripts/final.json'
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

    for i in range(len(xData)):
        xfinal.append(np.array(xData[i])[np.argsort(xData[i])])
        yfinal.append(np.array(yData[i])[np.argsort(xData[i])])

    plt.plot(xfinal[0], yfinal[0], '1', label='PER=' + '% 6.3f' % pers[0], color='b')
    #plt.plot(xfinal[0], np.power(xfinal[0],2), '1', label='PER=' + '% 6.3f' % pers[0], color='b')
    plt.plot(xfinal[1], yfinal[1], '2', label='PER=' + '% 6.3f' % pers[1], color='g')
    plt.plot(xfinal[2], yfinal[2], '3', label='PER=' + '% 6.3f' % pers[2], color='m')
    plt.plot(xfinal[3], yfinal[3], '4', label='PER=' + '% 6.3f' % pers[3], color='y')

    optimizedParameters, pcov = opt.curve_fit(LinFun, xData[0], yData[0])
    optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, xData[0], yData[0])
    plt.plot(xData[0], LinFun(xData[0], *optimizedParameters),linestyle='--', color='b',  label='O(n)')
    plt.plot(xData[0], QuadFun(xData[0], *optimizedParameters2),linestyle='--', color='r', label='O(n^2)')
    #
    # optimizedParameters, pcov = opt.curve_fit(LinFun, xData[1], yData[1])
    # optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, xData[1], yData[1])
    # plt.plot(xData[1], LinFun(xData[1], *optimizedParameters), color='b', label='O(n)')
    # plt.plot(xData[1], QuadFun(xData[1], *optimizedParameters2), color='r', label='O(n^2)')
    #
    # optimizedParameters, pcov = opt.curve_fit(LinFun, xData[2], yData[2])
    # optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, xData[2], yData[2])
    # plt.plot(xData[2], LinFun(xData[2], *optimizedParameters), color='b', label='O(n)')
    # plt.plot(xData[2], QuadFun(xData[2], *optimizedParameters2), color='r', label='O(n^2)')
    #
    # optimizedParameters, pcov = opt.curve_fit(LinFun, xData[3], yData[3])
    # optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, xData[3], yData[3])
    # plt.plot(xData[3], LinFun(xData[3], *optimizedParameters), color='b', label='O(n)')
    # plt.plot(xData[3], QuadFun(xData[3], *optimizedParameters2), color='r', label='O(n^2)')
    plt.legend()
    plt.grid()
    plt.show()
runtime()
