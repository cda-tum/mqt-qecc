import numpy as np
from bposd.hgp import hgp

from ldpc.codes import ring_code

# toric code
for i in range(2, 60):
    h1 = ring_code(i)
    h2 = ring_code(i)

    qcode = hgp(h1, h2, compute_distance=True)
    if qcode.test:  # to make sure it is a valid css code
        np.savetxt(f"examples/toricCodes/toric_{qcode.code_params}_hx.txt", qcode.hx, fmt='%d', newline='\n')
