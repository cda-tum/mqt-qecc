from mqt.qecc import *
import glob

code = Code("/home/luca/Documents/codeRepos/qecc/examples/hgp_(4,8)-[[5408,18,26]]_hx.txt", 18)
outpath = "./dp-sims-bindings.out"
outfile =  open(outpath, 'w')

runsPerEr = 1000000
currPer = 0.00001
maxPer = 0.1
stepSize = 0.00002
nrFailedRuns = 0
codeK = code.K

while(currPer < maxPer):
    nrFailedRuns = 0
    for i in range(runsPerEr):
        decoder = ImprovedUFD()
        decoder.set_code(code)
        err = sample_iid_pauli_err(code.N, currPer)
        decoder.decode(code.get_syndrome(err))
        result = decoder.result.estimate
        if(not code.is_stabilizer(result)):
            nrFailedRuns += 1
        currPer += stepSize
    fracFailed = nrFailedRuns/runsPerEr
    wer = fracFailed/codeK
    outfile.write(str(currPer) + ":" + str(wer) + '\n')
    currPer += stepSize