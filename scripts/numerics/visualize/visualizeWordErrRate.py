import json
import numpy as np

import matplotlib.pyplot as plt


def wer():
    plt.rcParams.update({'font.size': 14})
    inputFilename = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/decodingPerformance/1024_lp_code/dp-heur-singlesmallest-100k-runs.json'
    inputFilename2 = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/decodingPerformance/1024_lp_code/dp-original-stdgrowth-100k-runs.json'
    inputFilename3 = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/decodingPerformance/1024_lp_code/dp-heur-standardg-100k-runs.json'
    fig,ax = plt.subplots(1)

    f = open(inputFilename)
    data = json.load(f)
    xData = []
    yData = []
    for key in data:
        xData.append(float(key))
        yData.append(data[key])

    order = np.argsort(xData)
    xDataf =np.array(xData)[order]
    yDataf =np.array(yData)[order]

    f2 = open(inputFilename2)
    data2 = json.load(f2)
    xData2 = []
    yData2 = []
    for key2 in data2:
        xData2.append(float(key2))
        yData2.append(data2[key2])

    order2 = np.argsort(xData2)
    xData2f =np.array(xData2)[order2]
    yData2f =np.array(yData2)[order2]

    f3 = open(inputFilename3)
    data3 = json.load(f3)
    xData3 = []
    yData3 = []
    for key3 in data3:
        xData3.append(float(key3))
        yData3.append(data3[key3])

    order3 = np.argsort(xData3)
    xData3f =np.array(xData3)[order3]
    yData3f =np.array(yData3)[order3]

    ax.plot(xDataf, yDataf, '-x', label='heuristic SSG', color='b')
    ax.plot(xData2f, yData2f, '-o', label='general UFD', color='g')
    ax.plot(xData3f, yData3f, '-x', label='heuristic AG', color='r')
    print(yData)
    print(yData2)
    print(yData3)

    ax.set_xlabel('physical X-error rate')
    ax.set_ylabel('WER')
    ax.legend()
    ax.grid()
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.gca().set_position([0, 0, 1, 1])
    handles,labels = ax.get_legend_handles_labels()

    handles = [handles[1], handles[2], handles[0] ]
    labels = [labels[1], labels[2], labels[0] ]

    ax.legend(handles,labels,loc=4)
    plt.show()
def werComp():
    plt.rcParams.update({'font.size': 15})
    inputFilename = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/decodingPerformance/1024_lp_code/dp-heur-singlerandomg-100k-runs.json'
    inputFilename2 = '/home/luca/Documents/codeRepos/qecc/scripts/numerics/data/decodingPerformance/1024_lp_code/dp-original-stdgrowth-100k-runs.json'
    fig,ax = plt.subplots(1)

    f = open(inputFilename)
    data = json.load(f)
    xData = []
    yData = []
    for key in data:
        xData.append(float(key))
        yData.append(data[key])

    f2 = open(inputFilename2)
    data2 = json.load(f2)
    xData2 = []
    yData2 = []
    for key2 in data2:
        xData2.append(float(key2))
        yData2.append(data2[key2])

    ax.plot(xData, yData, '-d', label='heuristic', color='b')
    ax.plot(xData2, yData2, '-o', label='decoder', color='g')
    print(yData)
    print(yData2)

    ax.xlabel('physical X-error rate')
    ax.ylabel('WER')
    ax.legend()
    ax.grid()
    ax.xscale('log')
    ax.yscale('log')
    ax.gca().set_position([0, 0, 1, 1])
    ax.savefig("werComp.svg")
    handles,labels = ax.get_legend_handles_labels()

    handles = [handles[2], handles[1], handles[0]]
    labels = [labels[2], labels[1], labels[0]]

    ax.legend(handles,labels,loc=2)
    plt.show()

wer()
#werComp()
