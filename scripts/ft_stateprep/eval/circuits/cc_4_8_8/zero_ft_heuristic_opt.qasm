OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
qreg z_anc[5];
qreg x_anc[3];
qreg a12[3];
qreg a13[1];
qreg a14[1];
creg z_c[5];
creg x_c[3];
creg c12[3];
creg c13[1];
creg c14[1];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[9];
h q[10];
cx q[4],q[8];
cx q[10],q[15];
cx q[4],q[6];
cx q[9],q[15];
cx q[5],q[12];
cx q[4],q[7];
cx q[3],q[14];
cx q[2],q[6];
cx q[14],q[16];
cx q[10],q[13];
cx q[9],q[11];
cx q[8],q[5];
cx q[3],q[7];
cx q[1],q[4];
cx q[0],q[2];
cx q[12],q[9];
cx q[1],q[13];
cx q[0],q[3];
cx q[15],q[8];
cx q[6],q[16];
cx q[5],q[11];
cx q[4],q[10];
cx q[2],q[14];
cx q[1],z_anc[0];
cx q[5],z_anc[0];
cx q[9],z_anc[0];
cx q[10],z_anc[0];
cx q[15],z_anc[0];
cx q[0],z_anc[1];
cx q[2],z_anc[1];
cx q[7],z_anc[1];
cx q[16],z_anc[1];
cx q[0],z_anc[2];
cx q[2],z_anc[2];
cx q[7],z_anc[2];
cx q[16],z_anc[2];
cx q[1],z_anc[3];
cx q[3],z_anc[3];
cx q[4],z_anc[3];
cx q[6],z_anc[3];
cx q[14],z_anc[3];
cx q[1],z_anc[4];
cx q[4],z_anc[4];
cx q[9],z_anc[4];
cx q[11],z_anc[4];
measure z_anc[0] -> z_c[0];
measure z_anc[1] -> z_c[1];
measure z_anc[2] -> z_c[2];
measure z_anc[3] -> z_c[3];
measure z_anc[4] -> z_c[4];
h x_anc[0];
cx x_anc[0],q[0];
cx x_anc[0],a12[0];
cx x_anc[0],q[1];
cx x_anc[0],a12[1];
cx x_anc[0],q[8];
cx x_anc[0],q[10];
cx x_anc[0],a12[2];
cx x_anc[0],q[11];
cx x_anc[0],q[12];
cx x_anc[0],a12[0];
measure a12[0] -> c12[0];
cx x_anc[0],q[14];
cx x_anc[0],a12[2];
measure a12[2] -> c12[2];
cx x_anc[0],a12[1];
measure a12[1] -> c12[1];
cx x_anc[0],q[16];
h x_anc[0];
measure x_anc[0] -> x_c[0];
h x_anc[1];
cx x_anc[1],q[5];
cx x_anc[1],a13[0];
cx x_anc[1],q[10];
cx x_anc[1],q[12];
cx x_anc[1],a13[0];
measure a13[0] -> c13[0];
cx x_anc[1],q[13];
h x_anc[1];
measure x_anc[1] -> x_c[1];
h x_anc[2];
cx x_anc[2],q[1];
cx x_anc[2],a14[0];
cx x_anc[2],q[4];
cx x_anc[2],q[10];
cx x_anc[2],a14[0];
measure a14[0] -> c14[0];
cx x_anc[2],q[13];
h x_anc[2];
measure x_anc[2] -> x_c[2];
