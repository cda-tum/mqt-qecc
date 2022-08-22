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
    inputFilen = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/heuristic_growth-all_runs-100k.json'
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
            perXData.append(float(c))
            perYData.append(float(data[per][c]) / 1000.0)
        pers.append(float(per))
        xData.append(perXData)
        yData.append(perYData)

    xfinal = []
    yfinal = []
    colors = mcolors.BASE_COLORS

    for i in range(len(xData)):
        col, val = colors.popitem()
        if (col == 'w' or col == 'k'):
            col, val = colors.popitem()
            if (col == 'w' or col == 'k'):
                col, val = colors.popitem()
        label = '% 6.3f' % pers[i]
        orders.append(np.argsort(xData[i]))
        xfinal.append(np.array(xData[i])[orders[i]])
        yfinal.append( np.array(yData[i])[orders[i]])
        plt.plot(xfinal[i], yfinal[i], 'o', label='p=' + label, color=col)
        optimizedParameters, pcov = opt.curve_fit(LinFun, xfinal[i], yfinal[i])
        #optimizedParameters2, pcov2 = opt.curve_fit(QuadFun, xData[i][17:], yData[i][17:])
        #optimizedParameters2, pcov2 = opt.curve_fit(PowThreeFun, xData[i][17:], yData[i][17:])
        plt.plot(xfinal[i][10:], LinFun(xfinal[i][10:], *optimizedParameters),'--', color=col, label='O(n)')
        #plt.plot(xData[i][17:], PowThreeFun(xData[i][17:], *optimizedParameters2), '--',color=col, label='O(n^3)')
        #print(str(optimizedParameters[0]) + 'x+' + str(optimizedParameters[1]))
        print(xfinal)
        print(yfinal)


    plt.legend()
    plt.ylabel("avg runtime (s) to decode 100k samples")
    plt.xlabel("code length n")
    plt.grid()
    plt.show()

runtime()
