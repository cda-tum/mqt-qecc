from pathlib import Path

from mqt import qecc

code = qecc.Code("../../examples/lp_(4,8)-[[1024,18,nan]]_hx.txt", "../../examples/lp_(4,8)-[[1024,18,nan]]_hz.txt")
code.K = 18
outpath = "./dp-sims-bindings.out"

runsPerEr = 1
currPer = 0.00001
maxPer = 0.00003
stepSize = 0.00001
nrFailedRuns = 0
codeK = code.K

with Path(outpath).open("w") as outfile:
    while currPer < maxPer:
        nrFailedRuns = 0
        for _ in range(runsPerEr):
            decoder = qecc.UFHeuristic()
            decoder.set_code(code)
            err = qecc.sample_iid_pauli_err(code.N, currPer)
            decoder.decode(code.get_x_syndrome(err))
            result = decoder.result.estimate
            if not code.is_stabilizer(result):
                nrFailedRuns += 1
            currPer += stepSize
        fracFailed = nrFailedRuns / runsPerEr
        wer = fracFailed / codeK
        res = str(currPer) + ":" + str(wer) + "\n"
        print(res)
        outfile.write(res)
        currPer += stepSize
