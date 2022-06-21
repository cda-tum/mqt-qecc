import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd


def runtime():
    inputFilenames = ['/home/luca/Documents/uf-simulations/testrun/out09-06-2022.json']
    fig, ax = plt.subplots()
    colors = mcolors.BASE_COLORS
    xData = []
    yData = []

    for file in inputFilenames:
        data = pd.read_json(file)['runs']
        fileXData = []
        fileYData = []

        for r in range(0, len(data)):
            currRun = data[r]['run']
            physErrRate = currRun['physicalErrRate']
            fileXData.append(physErrRate)
            rData = currRun['data']
            avgDecodingTime = 0
            for decodingRun in rData:
                avgDecodingTime += decodingRun['decodingTime(ms)']
            if avgDecodingTime == 0:
                avgDecodingTime = 0
            else:
                avgDecodingTime = avgDecodingTime / len(rData)
            fileYData.append(avgDecodingTime)
        xData.append(fileXData)
        yData.append(fileYData)

    for i in range(0, len(xData)):
        col, val = colors.popitem()
        if (col == 'w' or col == 'k'):
            col, val = colors.popitem()
            if (col == 'w' or col == 'k'):
                col, val = colors.popitem()
        ax.plot(xData[i], yData[i], label='decoder ' + str(i), color=col)
    ax.set_xlabel('physical error rate')
    ax.set_ylabel('Avg runtime(ms)')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()


runtime()
