"""Constructs the lifted parity check codes for the lifted parity check codes."""
# using tools from https://github.com/quantumgizmos/bp_osd
import ldpc.protograph as pt
import numpy as np
from lifted_hgp import LiftedHgp

a1 = pt.array([[(0), (11), (7), (12)], [(1), (8), (1), (8)], [(11), (0), (4), (8)], [(6), (2), (4), (12)]])

qcode = LiftedHgp(lift_parameter=32, a=a1, b=a1)
if qcode.test:
    np.savetxt(f"./lp_{qcode.code_params}_hz.txt", qcode.hz, fmt="%d", newline="\n")
    np.savetxt(f"./lp_{qcode.code_params}_hx.txt", qcode.hx, fmt="%d", newline="\n")
