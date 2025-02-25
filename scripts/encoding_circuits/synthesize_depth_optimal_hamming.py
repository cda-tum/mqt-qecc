from __future__ import annotations

import numpy as np
from qiskit import qasm2

from mqt.qecc import CSSCode
from mqt.qecc.circuit_synthesis import depth_optimal_encoding_circuit

hamming = np.array([
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
    [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
])

logicals = np.array([
    [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],
    [1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
    [1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0],
    [1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
    [1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0],
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0],
])

hamming_code = CSSCode(3, hamming, hamming)
hamming_code.Lx = logicals
hamming_code.Lz = logicals

qc, enc = depth_optimal_encoding_circuit(hamming_code, min_depth=8, max_depth=12, min_timeout=300, max_timeout=3600)

qc_str = qasm2.dumps(qc)
with open("hamming_depth_optimal.qasm", "w", encoding="utf-8") as f:
    f.write(str(list(enc)))
    f.write(qc_str)
