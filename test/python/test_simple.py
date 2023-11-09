"""Test simple example with hamming code."""

from __future__ import annotations

import numpy as np

from mqt.qecc import Code, UFHeuristic, sample_iid_pauli_err


def test_basic() -> None:
    """Test basic functionality with Hamming code."""
    h = [
        [True, False, False, True, False, True, True],
        [False, True, False, True, True, False, True],
        [False, False, True, False, True, True, True],
    ]
    code = Code(h, h)
    decoder = UFHeuristic()
    decoder.set_code(code)
    x_err = sample_iid_pauli_err(code.n, 0.05)
    decoder.decode(code.get_x_syndrome(x_err))
    result = decoder.result
    residual_err = np.array(x_err) ^ np.array(result.estimate)

    print(result)
    print(code.is_x_stabilizer(residual_err))
    print(np.array(x_err).astype(int))
