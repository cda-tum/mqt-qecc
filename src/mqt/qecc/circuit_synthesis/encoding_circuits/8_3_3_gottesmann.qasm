// Messaging qubits: 1,2,4
// Stabilizers: ["XXXXXXXX", "ZZZZZZZZ", "IXIXYZYZ", "IXZYIXZY", "IYXZXZIY"]
// x_logicals=["XXIIIZIZ", "XIXZIIZI", "XIIZXZII"]
// z_logicals=["IZIZIZIZ", "IIZZIIZZ", "IIIIZZZZ"]
// Note that the encoding circuit does not preserve the sign of the stabilizers
// This needs to be accounted for when using this circuit for encoding/decoding
OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
h q[3];
sx q[7];
h q[6];
sx q[5];
cz q[3],q[4];
s q[2];
s q[1];
sx q[0];
cx q[7],q[4];
cx q[6],q[2];
cx q[5],q[3];
cz q[0],q[1];
cx q[7],q[1];
s q[6];
s q[2];
cx q[0],q[5];
cx q[5],q[7];
cx q[3],q[6];
s q[1];
cz q[0],q[2];
cx q[7],q[2];
cx q[6],q[4];
s q[5];
cx q[3],q[1];
h q[0];

z q[3];
z q[7];
x q[0];
