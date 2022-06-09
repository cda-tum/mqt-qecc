import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import json
import numpy as np


def wer():
    inputFilename = '/home/luca/Documents/uf-simulations/testrun/data09-06-2022.json'

    fig, ax = plt.subplots()

    f = open(inputFilename)
    data = json.load(f)
    xData = []
    yData = []

    for key in data:
        xData.append(float(key))
        yData.append(data[key])

    n = 10
    ax.plot(xData[::n], yData[::n], 'o', label='code 1', color='b')
    ax.set_xlabel('physical X-error rate')
    ax.set_ylabel('WER')
    ax.legend()

    ax.set_xscale('log')
    ax.set_yscale('log')
    print(xData)
    print(yData)
    plt.show()

wer()