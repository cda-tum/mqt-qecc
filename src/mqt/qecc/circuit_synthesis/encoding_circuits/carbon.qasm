# Encoding Qubits: 0,4
OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
h q[5];
h q[8];
h q[9];
h q[10];
h q[11];
cx q[9],q[2];
cx q[5],q[2];
cx q[11],q[1];
cx q[8],q[7];
cx q[5],q[4];
cx q[10],q[3];
cx q[9],q[1];
cx q[2],q[0];
cx q[11],q[7];
cx q[0],q[6];
cx q[1],q[4];
cx q[8],q[3];
cx q[4],q[7];
cx q[10],q[6];
cx q[3],q[2];
cx q[1],q[0];
cx q[6],q[11];
cx q[5],q[10];
cx q[7],q[9];
cx q[0],q[8];
cx q[2],q[1];
cx q[4],q[3];
