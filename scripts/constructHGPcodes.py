import numpy as np
from bposd.hgp import hgp

from ldpc.codes import ring_code

# med sized HGP code from

h = np.loadtxt("examples/mkmn_24_6_10.txt").astype(int)
qcode = hgp(h)  # construct quantum LDPC code using the symmetric hypergraph product
seed_code = np.loadtxt(f"examples/mkmn_24_6_10.txt").astype(int)
# print(seed_code)
qcode = hgp(seed_code, compute_distance=True)
qcode.canonical_logicals()
qcode.test()
print(qcode.code_params)
print(qcode.hx)
print("hz:")
print(qcode.hz)
np.savetxt(f"examples/hgp_{qcode.code_params}_hx.txt", qcode.hx, fmt='%d', newline='\n')

# toric code
for i in range(2, 1000):
    h1 = ring_code(i)
    h2 = ring_code(i)

    qcode = hgp(h1, h2, compute_distance=True)
    if qcode.test:  # to make sure it is a valid css code
        np.savetxt(f"examples/toric_{qcode.code_params}_hx.txt", qcode.hx, fmt='%d', newline='\n')
