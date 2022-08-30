from mqt.qecc import *

codes = [
    "../../examples/toricCodes/toric-[[128,2,8]]_hx.txt",
    "../../examples/toricCodes/toric-[[162,2,9]]_hx.txt",
    "../../examples/toricCodes/toric-[[242,2,11]]_hx.txt",
    "../../examples/toricCodes/toric-[[338,2,13]]_hx.txt",
    "../../examples/toricCodes/toric-[[392,2,14]]_hx.txt",
    "../../examples/toricCodes/toric-[[450,2,15]]_hx.txt",
    "../../examples/toricCodes/toric-[[512,2,16]]_hx.txt",
    "../../examples/toricCodes/toric-[[578,2,17]]_hx.txt",
    "../../examples/toricCodes/toric-[[648,2,18]]_hx.txt",
    "../../examples/toricCodes/toric-[[722,2,19]]_hx.txt",
    "../../examples/toricCodes/toric-[[882,2,21]]_hx.txt",
    "../../examples/toricCodes/toric-[[968,2,22]]_hx.txt",
    "../../examples/toricCodes/toric-[[1058,2,23]]_hx.txt"
]
outpath = "./rt-sims-bindings.out"
outfile =  open(outpath, 'w')

nrSamples = 50
nrRuns = 1000
per = 0.01

for codePath in codes:
    sampleSum = 0.0
    for i in range(nrSamples):
        runsSum = 0.0
        for j in range(nrRuns):
            code = Code(codePath)
            err = sample_iid_pauli_err(code.N, per)
            decoder = UFHeuristic()
            decoder.set_code(code)
            syndr = code.get_syndrome(err)
            decoder.decode(syndr)
            time = decoder.result.decoding_time
            runsSum += time
        sampleSum += runsSum
    outp = codePath + ":" + str((sampleSum / nrSamples))
    print(outp)
    outfile.write(outp + '\n')
    outfile.flush()