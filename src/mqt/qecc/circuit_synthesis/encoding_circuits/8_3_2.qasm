// Messaging qubits: 0, 1, 3

OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
h q[7];
cx q[7],q[3];
cx q[1],q[5];
cx q[0],q[4];
cx q[7],q[5];
cx q[4],q[6];
cx q[0],q[2];
cx q[3],q[1];
cx q[7],q[6];
cx q[5],q[4];
cx q[3],q[2];
cx q[1],q[0];
