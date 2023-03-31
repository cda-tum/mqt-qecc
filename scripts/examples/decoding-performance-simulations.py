"""Simulate the decoding performance of the UF heuristic decoder."""
from __future__ import annotations

from pathlib import Path

from mqt import qecc

code = qecc.Code("../../examples/lp_(4,8)-[[1024,18,nan]]_hx.txt", "../../examples/lp_(4,8)-[[1024,18,nan]]_hz.txt")
code.K = 18
outpath: str = "./dp-sims-bindings.out"

runs_per_p: int = 1
curr_p: float = 0.00001
max_per: float = 0.00003
step_size: float = 0.00001
nr_failed_runs: int = 0
code_k: int = code.K

with Path(outpath).open("w") as outfile:
    while curr_p < max_per:
        nr_failed_runs = 0
        for _ in range(runs_per_p):
            decoder = qecc.UFHeuristic()
            decoder.set_code(code)
            err = qecc.sample_iid_pauli_err(code.N, curr_p)
            decoder.decode(code.get_x_syndrome(err))
            result = decoder.result.estimate
            if not code.is_stabilizer(result):
                nr_failed_runs += 1
        frac_failes = nr_failed_runs / runs_per_p
        wer = frac_failes / code_k
        res = str(curr_p) + ":" + str(wer) + "\n"
        print(res)
        outfile.write(res)
        curr_p += step_size
