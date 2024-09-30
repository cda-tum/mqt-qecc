from __future__ import annotations

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
import typing

from pathlib import Path
import os
import numpy as np
import pandas as pd

from qiskit import QuantumCircuit
from mqt.qecc import codes as codes
import mqt.qecc.ft_stateprep as ftsp

import qsample as qs

from concurrent.futures import ProcessPoolExecutor
from time import time

if typing.TYPE_CHECKING:
    import numpy.typing as npt

# The directory where the circuits are stored
# The expected directory structure is:
# circ_dir
# ├── code_name
# │   ├── {zero,plus}_{heuristic,opt}.qasm
prefix = (Path(__file__) / "../circuits/").resolve()
# Synthesis parameters
max_timeout = 3000
min_timeout = 8
max_ancillas = 8

# Simulation parameters
err_params = {'q': [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]}
err_model = qs.noise.E1_1
shots_dss = 8000
p_max = {"q" : 0.01}
L = 3

# Helper functions
def _extract_parameters(file_path):
    # Extract the code name from the parent directory
    code_name = os.path.basename(os.path.dirname(file_path))
    zero_state = 'zero' in file_path
    # Extract the procedure from the filename
    if 'heuristic' in file_path:
        procedure = 'heuristic'
    elif 'opt' in file_path:
        procedure = 'opt'
    else:
        procedure = 'unknown'
    
    return code_name, zero_state, procedure

def _codes_from_matrix(code_name):
    if code_name == '11_1_3':
        matrix = np.array([[1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
                 [0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                 [0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0],
                 [0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0],
                 [0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1]])
        return codes.CSSCode(distance=3, Hx=matrix, Hz=matrix)
    elif code_name == '16_2_4':
        matrix = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1],
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1],
                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1],
                 [0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0],
                 [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                 [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0]])
        return codes.CSSCode(distance=4, Hx=matrix, Hz=matrix)
    elif code_name == "hypercube":
        matrix = np.array([[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
                        [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
                        [1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1],
                        [0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0],
                        [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0]])
        return codes.CSSCode(distance=4, Hx=matrix, Hz=matrix)



def run_parallel(path: str):
    # load code and circuit
    code_name, zero_state, procedure = _extract_parameters(path)
    print(code_name, zero_state, procedure)
    if code_name in ['11_1_3', '16_2_4', 'hypercube']:
        code = _codes_from_matrix(code_name)
    else:
        code = codes.CSSCode.from_code_name(code_name)
    state_prep_circ = ftsp.StatePrepCircuit(QuantumCircuit.from_qasm_file(path), code, zero_state)

    starttime = time()

    verify_helper = ftsp.DeterministicVerificationHelper(state_prep_circ)
    verify_heuristic = verify_helper.get_solution(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancillas, use_optimal_verification=False)
    verify_optimal = verify_helper.get_solution(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancillas, use_optimal_verification=True)
    verify_global = verify_helper.get_global_solution(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancillas, max_possibilities=2)


    sim_heuristic = ftsp.NoisyDFTStatePrepSimulator(state_prep_circ.circ, verify_heuristic, code, err_model, zero_state)
    stats_heuristic = sim_heuristic.dss_logical_error_rates(err_params=err_params, p_max=p_max, L=L, shots=shots_dss)
    sim_optimal = ftsp.NoisyDFTStatePrepSimulator(state_prep_circ.circ, verify_optimal, code, err_model, zero_state)
    stats_optimal = sim_optimal.dss_logical_error_rates(err_params=err_params, p_max=p_max, L=L, shots=shots_dss)
    sim_global = ftsp.NoisyDFTStatePrepSimulator(state_prep_circ.circ, verify_global, code, err_model, zero_state)
    stats_global = sim_global.dss_logical_error_rates(err_params=err_params, p_max=p_max, L=L, shots=shots_dss)

    results = pd.DataFrame()
    verifications = ["heuristic", "optimal", "global"]
    for name, verify, stats in zip(verifications, [verify_heuristic, verify_optimal, verify_global], [stats_heuristic, stats_optimal, stats_global]):
        results = pd.concat([results, pd.DataFrame({'code': [code_name], 'zero_state': [zero_state], 'global_opt': [False], 'verification': [name], 'procedure': [procedure], 'verification_stabs_0': [verify[0].stabs], 'recovery_stabs_0': [verify[0].det_correction], 'flags_0': [verify[0].hook_corrections], 'verification_stabs_1': [verify[1].stabs], 'recovery_stabs_1': [verify[1].det_correction], 'flags_1': [verify[1].hook_corrections], 'logical_error_rates': [stats], 'time': [time()-starttime]})])
    
    results.to_csv(f"results_{code_name}_{zero_state}_{procedure}.csv")

    return results




if __name__ == '__main__':
    # get all filepaths in the directory
    file_paths = []
    for root, dirs, files in os.walk(prefix):
        for file in files:
            if file.endswith(".qasm"):
                file_paths.append(os.path.join(root, file))
    
    # run in parallel
    with ProcessPoolExecutor() as executor:
        dfs = list(executor.map(run_parallel, file_paths))

    results = pd.concat(dfs, ignore_index=True)
    results.to_csv('results.csv')
    

