// messaging qubits: 4

OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
sx q[0];
h q[1];
h q[2];
cx q[2],q[4];
sx q[2];
cz q[1],q[2];
s q[1];
sx q[2];
cx q[4],q[3];
sx q[3];
cz q[0],q[3];
h q[0];
cx q[0],q[2];
sx q[3];
cx q[3],q[1];
z q[4];
x q[3];
z q[2];
y q[0];
y q[1];
