OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
qreg x_anc[1];
creg x_c[1];
h q[0];
h q[3];
h q[4];
h q[5];
cx q[3],q[1];
cx q[4],q[3];
cx q[4],q[6];
cx q[0],q[2];
cx q[5],q[0];
cx q[5],q[6];
cx q[2],q[3];
cx q[3],q[5];
h x_anc[0];
cx x_anc[0],q[2];
cx x_anc[0],q[4];
cx x_anc[0],q[5];
h x_anc[0];
measure x_anc[0] -> x_c[0];