from pathlib import Path

from mqt import qecc

codes = [
    "../../examples/toricCodes/toric_(nan,nan)-[[8,2,2]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[18,2,3]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[32,2,4]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[50,2,5]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[72,2,6]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[98,2,7]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[128,2,8]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[162,2,9]]_hz.txt",
    "../../examples/toricCodes/toric_(nan,nan)-[[200,2,10]]_hz.txt",
]
outpath = "./rt-sims-bindings.out"

nrSamples = 1
nrRuns = 1
per = 0.01

with Path(outpath).open("w") as outfile:
    for codePath in codes:
        sampleSum = 0.0
        for _ in range(nrSamples):
            runsSum = 0.0
            for _ in range(nrRuns):
                code = qecc.Code(codePath)
                err = qecc.sample_iid_pauli_err(code.N, per)
                decoder = qecc.UFHeuristic()
                decoder.set_code(code)
                syndr = code.get_x_syndrome(err)
                decoder.decode(syndr)
                time = decoder.result.decoding_time
                runsSum += time
            sampleSum += runsSum
        outp = codePath + ":" + str(sampleSum / nrSamples)
        print(outp)
        outfile.write(outp + "\n")
        outfile.flush()
