import json

import matplotlib.pyplot as plt


def wer():
    plt.rcParams.update({'font.size': 15})
    inputFilename = '/home/luca/Documents/codeRepos/qecc/scripts/subm/dp-10k-0001-01.json'
    inputFilename2 = '/home/luca/Documents/codeRepos/qecc/scripts/subm/dp-10k-0001-01-singlecg.json'
    #inputFilename = '/home/luca/Documents/codeRepos/qecc/scripts/subm/dp-impr-10k-0001-01.json'
    #inputFilename2 = '/home/luca/Documents/codeRepos/qecc/scripts/subm/dp-impr-sscg.json'

    fig, ax = plt.subplots()

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

    ax.plot(xData, yData, 'd', label='standard growth', color='b')
    ax.plot(xData2, yData2, 'o', label='single smallest growth', color='r')
    ax.set_xlabel('physical X-error rate')
    ax.set_ylabel('WER')
    ax.legend()
    plt.grid()
    ax.set_xscale('log')

    print(xData)
    print(yData)
    plt.show()
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("test.svg")
wer()
