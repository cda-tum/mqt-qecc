import ldpc.protograph as pt
import numpy as np

from lifted_hgp import lifted_hgp

a1 = pt.array([
    [(0), (11), (7), (12)],
    [(1), (8), (1), (8)],
    [(11), (0), (4), (8)],
    [(6), (2), (4), (12)]])

qcode = lifted_hgp(lift_parameter=13, a=a1, b=a1)
if qcode.test:
    np.savetxt(f"examples/lp_{qcode.code_params}_hx.txt", qcode.hx, fmt='%d', newline='\n')
