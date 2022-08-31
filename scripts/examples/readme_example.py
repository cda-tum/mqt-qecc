from mqt.qecc import *
import numpy as np

code = Code('/path/to/Hx', 'path/to/Hz')
decoder = UFHeuristic()
decoder.set_code(code)
x_err = sample_iid_pauli_err(code.N, 0.05)
decoder.decode(code.get_x_syndrome(x_err))
result = decoder.result
print(result)
residual_err = np.array(x_err)^np.array(result.estimate)
print(code.is_x_stabilizer(residual_err))