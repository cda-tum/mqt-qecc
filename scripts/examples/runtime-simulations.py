from mqt.qecc import *

codes = [
   "../../examples/toricCodes/toric_(nan,nan)-[[8,2,2]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[18,2,3]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[32,2,4]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[50,2,5]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[72,2,6]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[98,2,7]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[128,2,8]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[162,2,9]]_hz.txt",
   "../../examples/toricCodes/toric_(nan,nan)-[[200,2,10]]_hz.txt"
]
outpath = "./rt-sims-bindings.out"
outfile =  open(outpath, 'w')

nrSamples = 1
nrRuns = 1
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
            syndr = code.get_x_syndrome(err)
            decoder.decode(syndr)
            time = decoder.result.decoding_time
            runsSum += time
        sampleSum += runsSum
    outp = codePath + ":" + str((sampleSum / nrSamples))
    print(outp)
    outfile.write(outp + '\n')
    outfile.flush()