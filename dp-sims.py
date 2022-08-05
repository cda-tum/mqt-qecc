import os
import json
import glob
import os

nrRuns = 2
nrSamples = 2
per = 0.001
#outpath = "/home/luca/Documents/uf-simulations/runtime/original/"
outpath = "/home/berent/ufpaper/simulations/montecarlo/fixed/out/"

for cnt in range(nrSamples):
    for i in range(nrRuns):
        os.system("parallel ./qecc_app " + str(per) + " " + outpath + "< codenames.txt >>out"
                  + str(per) + ".log")

# Parse outputs
filepath = outpath
files = glob.glob(filepath+"*.txt")
outfilename = filepath+"parsed.json"
outfile =  open(outfilename, 'w')
outfile.write("{ \"" + str(per) + "\":")
codeData = {}

for file in files:
    with open(file, 'r') as original: data = original.readlines()
    code = file.split(".")[0].split("/")
    code = code[len(code)-1]

    cnt = 0
    sumRun = 0
    samples = []
    for line in data:
        if cnt < nrRuns:
            sumRun += float(line)
            cnt = cnt + 1
        else:
            samples.append(sumRun)
            sumRun = 0
            cnt = 0

    sum = 0
    for t in samples:
        sum += t
    codeData[code] = (sum/nrSamples)

outfile.write(json.dumps(codeData)+",")
outfile.close()

with open(outfilename, 'rb+') as filehandle:
    filehandle.seek(-1, os.SEEK_END)
    filehandle.truncate()

outfile =  open(outfilename, 'a')
outfile.write("}")
