from mqt.qecc import *
import glob

codes = glob.glob("/home/luca/Documents/codeRepos/qecc/examples/toricCodes2/*.txt")
outpath = "./rt-sims-bindings.out"
outfile =  open(outpath, 'w')

nrSamples = 100000
nrRuns = 50
per = 0.01

for codePath in codes:
    sampleSum = 0.0
    for i in range(nrSamples):
        runsSum = 0.0
        for j in range(nrRuns):
            code = Code(codePath)
            err = sample_iid_pauli_err(code.N, per)
            decoder = ImprovedUFD()
            decoder.set_code(code)
            syndr = code.get_syndrome(err)
            decoder.decode(syndr)
            time = decoder.result.decoding_time
            runsSum += time
            decoder.reset()
        sampleSum += runsSum
    outp = codePath + ":" + str((sampleSum / nrSamples))
    print(outp)
    outfile.write(outp + '\n')
    outfile.flush()