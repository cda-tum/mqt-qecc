from __future__ import annotations
from mqt.qecc import *
import numpy as np

def test_basic():
    H = [
        [1 ,0 ,0 ,1 ,0 ,1 ,1],
        [0 ,1 ,0 ,1 ,1 ,0 ,1],
        [0 ,0 ,1 ,0 ,1 ,1 ,1],
    ]
    code = Code(H, H)
    decoder = UFHeuristic()
    decoder.set_code(code)
    x_err = sample_iid_pauli_err(code.N, 0.05)
    decoder.decode(code.get_x_syndrome(x_err))
    result = decoder.result
    print(result)
    residual_err = np.array(x_err)^np.array(result.estimate)
    print(code.is_x_stabilizer(residual_err))