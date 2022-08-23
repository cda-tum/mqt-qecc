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
    inputFilen = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/runtime/rt-impr-100k-stdgrowth-02.json'
    inputFilen2 = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/runtime/rt-original-100k-stdgrowth-02.json'
    colors = mcolors.BASE_COLORS
    xData = []
    yData = []
    pers = []
    orders = []

    xData2 = []
    yData2 = []
    pers2 = []
    orders2 = []

    with open(inputFilen) as data_file:
        data = json.load(data_file)

    with open(inputFilen2) as data_file2:
        data2 = json.load(data_file2)

    for per in data:
        perXData = []
        perYData = []

        for c in data[per]:
            perXData.append(float(c))
            perYData.append(float(data[per][c]) / 1000.0)
        pers.append(float(per))
        xData.append(perXData)
        yData.append(perYData)

    for per in data2:
        perXData = []
        perYData = []

        for c in data2[per]:
            perXData.append(float(c))
            perYData.append(float(data2[per][c]) / 1000.0)
        pers2.append(float(per))
        xData2.append(perXData)
        yData2.append(perYData)

    xfinal = []
    yfinal = []
    cols = []
    colors = mcolors.BASE_COLORS
    colors2 = mcolors.BASE_COLORS

    for i in range(len(xData)):
        col, val = colors.popitem()
        if (col == 'w' or col == 'k'):
            col, val = colors.popitem()
            if (col == 'w' or col == 'k'):
                col, val = colors.popitem()
                cols.append(col)
        label = '% 6.3f' % pers[i]
        orders.append(np.argsort(xData[i]))
        xfinal.append(np.array(xData[i])[orders[i]])
        yfinal.append(np.array(yData[i])[orders[i]])
        plt.plot(xfinal[i], yfinal[i], 'o', label='p=' + label, color=col)
        print(xfinal)
        print(yfinal)

    start = 24
    optimizedParameters, pcov = opt.curve_fit(LinFun, xfinal[0][start:], yfinal[0][start:])
    plt.plot(xfinal[0][start:], LinFun(xfinal[0][start:], *optimizedParameters),'--', color=cols[0], label='O(n)')

    start2 = 7
    optimizedParameters2, pcov = opt.curve_fit(LinFun, xfinal[1][start2:], yfinal[1][start2:])
    plt.plot(xfinal[1][start2:], LinFun(xfinal[1][start2:], *optimizedParameters2),'--', color='m', label='O(n)')

    xfinal2 = []
    yfinal2 = []
    cols2=[]
    for i in range(len(xData2)):
        col, val = colors.popitem()
        cols2.append(col)
        if (col == 'w' or col == 'k'):
            col, val = colors.popitem()
            cols2.append(col)
            if (col == 'w' or col == 'k'):
                col, val = colors.popitem()
                cols2.append(col)
        label = '% 6.3f' % pers2[i]
        orders2.append(np.argsort(xData2[i]))
        xfinal2.append(np.array(xData2[i])[orders2[i]])
        yfinal2.append(np.array(yData2[i])[orders2[i]])
        plt.plot(xfinal2[i], yfinal2[i], 'o', label='p=' + label, color=col)

        #optimizedParameters2, pcov2 = opt.curve_fit(PowThreeFun, xfinal2[i], yfinal2[i])
        #plt.plot(xfinal2[i], PowThreeFun(xfinal2[i], *optimizedParameters2), '--',color=col, label='O(n^3)')
        print(xfinal2)
        print(yfinal2)

    #optimizedParameters3, pcov3 = opt.curve_fit(QuadFun, xfinal2[0], yfinal2[0])
    #plt.plot(xfinal2[0], QuadFun(xfinal2[0], *optimizedParameters3), '--',color=cols2[0], label='O(n^3)')

    plt.legend()
    plt.ylabel("avg runtime (s) to decode 100k samples")
    plt.xlabel("code length n")
    plt.grid()
    plt.show()

runtime()
