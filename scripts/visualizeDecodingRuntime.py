import json
import math

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


def LinFun(x, a, b):
    return np.multiply(a, x) + b


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
        ax.plot(xData[i], yData[i], 'o', label='PER' + label, color=col)
        optimizedParameters, pcov = opt.curve_fit(LinFun, xData[i], yData[i])
        func = str(math.ceil(optimizedParameters[1])) + 'x^2+' + str(math.ceil(optimizedParameters[0]))
        ax.plot(xData[i], LinFun(xData[i], *optimizedParameters), color='r')
        print(xData[i])
        print(yData[i])
    # fits

    ax.set_xlabel('code size')
    ax.set_ylabel('Avg runtime(ms)')
    ax.legend()
    plt.show()


runtime()
