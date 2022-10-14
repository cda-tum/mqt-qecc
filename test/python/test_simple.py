from __future__ import annotations

import numpy as np

from mqt.qecc import Code, UFHeuristic, sample_iid_pauli_err


def test_basic() -> None:
    H = [
        [1, 0, 0, 1, 0, 1, 1],
        [0, 1, 0, 1, 1, 0, 1],
        [0, 0, 1, 0, 1, 1, 1],
    ]
    code = Code(H, H)
    decoder = UFHeuristic()
    decoder.set_code(code)
    x_err = sample_iid_pauli_err(code.n, 0.05)
    decoder.decode(code.get_x_syndrome(x_err))
    result = decoder.result
    residual_err = np.array(x_err) ^ np.array(result.estimate)

    print(result)
    print(code.is_x_stabilizer(residual_err))
    print(np.array(x_err).astype(int))