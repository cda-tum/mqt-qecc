# Encoding Qubits: 0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
h q[3];
h q[5];
h q[6];
h q[7];
cx q[7],q[4];
cx q[3],q[0];
cx q[3],q[4];
cx q[5],q[2];
cx q[0],q[1];
cx q[7],q[8];
cx q[4],q[5];
cx q[1],q[2];
cx q[6],q[3];
