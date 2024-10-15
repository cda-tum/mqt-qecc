# Encoding qubits: 0, 1

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
h q[2];
cx q[1],q[3];
cx q[2],q[0];
cx q[2],q[3];
cx q[0],q[1];
