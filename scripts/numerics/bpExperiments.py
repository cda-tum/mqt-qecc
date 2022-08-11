import numpy as np
from bposd.css_decode_sim import css_decode_sim
from bposd.hgp import hgp

h = np.loadtxt("examples/codes/classical_seed_codes/mkmn_16_4_6.txt").astype(int)
qcode = hgp(h)  # construct quantum LDPC code using the symmetric hypergraph product

osd_options = {
    'error_rate': 0.05,
    'target_runs': 1000,
    'xyz_error_bias': [0, 0, 1],
    'output_file': 'test.json',
    'bp_method': "ms",
    'ms_scaling_factor': 0,
    'osd_method': "osd_cs2",
    'osd_order': 42,
    'channel_update': None,
    'seed': 42,
    'max_iter': 0,
    'output_file': "test.json"
}

lk = css_decode_sim(hx=qcode.hx, hz=qcode.hz, **osd_options)
