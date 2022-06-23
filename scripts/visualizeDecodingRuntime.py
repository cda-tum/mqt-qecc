import json

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


def runtime():
    inputFilen = '/home/luca/Documents/uf-simulations/testrun/raw-final23-06-2022.json'
    fig, ax = plt.subplots()
    colors = mcolors.BASE_COLORS
    xData = []
    yData = []
    pers = []

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
        ax.plot(xData[i], yData[i], label='PER' + label, color=col)

    ax.set_xlabel('code size')
    ax.set_ylabel('Avg runtime(ms)')
    ax.legend()
    plt.show()


runtime()
