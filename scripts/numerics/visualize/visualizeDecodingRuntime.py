import json
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'


def LinFun(x, a, b):
    return np.multiply(a, x) + b


def QuadFun(x, a, b):
    return np.multiply(a, np.power(x, 2)) + b


def PowThreeFun(x, a, b):
    return np.multiply(a, np.power(x, 3)) + b

def PowFourFun(x, a, b):
    return np.multiply(a, np.power(x, 4)) + b


def runtime():
    plt.rcParams.update({'font.size': 15})
    inputFilen = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/runtime/rt-original-1k-01.json'
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


def runtimeComparison():
    plt.rcParams.update({'font.size': 14})
    inputFilen = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/runtime/rt-heur-1k-50s.json'
    inputFilen2 = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/runtime/rt-original-1k-01.json'
    colors = mcolors.BASE_COLORS
    xData = []
    yData = []
    pers = []
    xData2 = []
    yData2 = []
    pers2 = 0.0
    orders = []
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

        for c in data2[per]:
            xData2.append(float(c))
            yData2.append(float(data2[per][c]) / 1000.0)
        pers2 = float(per)

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
        label = '%2.2f' % pers[i]
        orders.append(np.argsort(xData[i]))
        xfinal.append(np.array(xData[i])[orders[i]])
        yfinal.append(np.array(yData[i])[orders[i]])
        plt.plot(xfinal[i], yfinal[i], 'o', label='UFH, p=' + label, color=col)
        start = 0
        optimizedParameters, pcov = opt.curve_fit(LinFun, xfinal[i][start:], yfinal[i][start:])
        plt.plot(xfinal[i][start:], LinFun(xfinal[i][start:], *optimizedParameters), '--', color=col,
                 label=r'$O(n)$')
        print(xfinal)
        print(yfinal)

    # general qlpd decoder data
    label = '%2.2f' % pers2
    orders2 = np.argsort(xData2)
    xfinal2 = np.array(xData2)[orders2]
    yfinal2 = np.array(yData2)[orders2]
    plt.plot(xfinal2, yfinal2, 'd', label='GD, p=' + label, color='green')
    start = 0
    optimizedParameters2, pcov = opt.curve_fit(QuadFun, xfinal2[start:], yfinal2[start:])
    plt.plot(xfinal2[start:], QuadFun(xfinal2[start:], *optimizedParameters2), '--', color='green', label=r'$O(n^{2})$')


    plt.legend()
    plt.ylabel("avg runtime (s) to decode $10^3$ samples")
    plt.xlabel("code length n")
    plt.grid()
    plt.show()


runtimeComparison()
