from mqt.qecc import *

code = Code("/home/luca/Documents/codeRepos/qecc/examples/toricCodes/toric_(nan,nan)-[[1058,2,23]]_hx.txt")
err = sample_iid_pauli_err(1058, 0.1)
syndrome = code.get_syndrome(err)
decoder = ImprovedUFD()
decoder.set_code(code)

decoder.decode(syndrome)
print(decoder.result)
print(code.is_stabilizer(decoder.result.estimate))
