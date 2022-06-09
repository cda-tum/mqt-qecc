import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

import json
import numpy as np


def runtime():
    inputFilename = '/home/luca/Documents/uf-simulations/testrun/out09-06-2022.json'
    fig, ax = plt.subplots()

    data = pd.read_json(inputFilename)['runs']

    plt.show()

runtime()