OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
qreg z_anc[5];
qreg x_anc[4];
qreg a8[4];
qreg a9[2];
qreg a10[1];
qreg a11[2];
creg z_c[5];
creg x_c[4];
creg c8[4];
creg c9[2];
creg c10[1];
creg c11[2];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[10];
h q[12];
cx q[6],q[16];
cx q[4],q[13];
cx q[3],q[7];
cx q[12],q[17];
cx q[10],q[18];
cx q[7],q[11];
cx q[6],q[3];
cx q[4],q[8];
cx q[2],q[16];
cx q[1],q[15];
cx q[0],q[14];
cx q[11],q[12];
cx q[10],q[13];
cx q[8],q[9];
cx q[6],q[0];
cx q[5],q[2];
cx q[4],q[3];
cx q[1],q[7];
cx q[17],q[18];
cx q[14],q[15];
cx q[5],q[9];
cx q[16],q[4];
cx q[13],q[11];
cx q[12],q[10];
cx q[7],q[6];
cx q[2],q[8];
cx q[0],q[1];
cx q[0],z_anc[0];
cx q[2],z_anc[0];
cx q[7],z_anc[0];
cx q[9],z_anc[0];
cx q[12],z_anc[0];
cx q[13],z_anc[0];
cx q[15],z_anc[0];
cx q[16],z_anc[0];
cx q[18],z_anc[0];
cx q[4],z_anc[1];
cx q[5],z_anc[1];
cx q[9],z_anc[1];
cx q[16],z_anc[1];
cx q[0],z_anc[2];
cx q[1],z_anc[2];
cx q[3],z_anc[2];
cx q[5],z_anc[2];
cx q[7],z_anc[2];
cx q[8],z_anc[2];
cx q[16],z_anc[2];
cx q[1],z_anc[3];
cx q[6],z_anc[3];
cx q[11],z_anc[3];
cx q[13],z_anc[3];
cx q[14],z_anc[3];
cx q[1],z_anc[4];
cx q[2],z_anc[4];
cx q[4],z_anc[4];
cx q[9],z_anc[4];
cx q[15],z_anc[4];
measure z_anc[0] -> z_c[0];
measure z_anc[1] -> z_c[1];
measure z_anc[2] -> z_c[2];
measure z_anc[3] -> z_c[3];
measure z_anc[4] -> z_c[4];

h x_anc[0];

cx x_anc[0],q[0];

cx x_anc[0],a8[0];

cx x_anc[0],q[2];

cx x_anc[0],a8[1];

cx x_anc[0],q[3];

cx x_anc[0],q[4];

cx x_anc[0],a8[2];

cx x_anc[0],q[6];

cx x_anc[0],q[9];

cx x_anc[0],a8[3];

cx x_anc[0],q[11];

cx x_anc[0],a8[2];
measure a8[2] -> c8[2];

cx x_anc[0],q[12];

cx x_anc[0],a8[0];
measure a8[0] -> c8[0];

cx x_anc[0],q[15];

cx x_anc[0],a8[3];
measure a8[3] -> c8[3];

cx x_anc[0],a8[1];
measure a8[1] -> c8[1];

cx x_anc[0],q[18];

h x_anc[0];
measure x_anc[0] -> x_c[0];

h x_anc[1];
cx x_anc[1],q[0];
cx x_anc[1],a9[0];
cx x_anc[1],q[2];
cx x_anc[1],a9[1];
cx x_anc[1],q[3];
cx x_anc[1],q[7];
cx x_anc[1],a9[0];
measure a9[0] -> c9[0];
cx x_anc[1],q[8];
cx x_anc[1],a9[1];
measure a9[1] -> c9[1];
cx x_anc[1],q[15];
h x_anc[1];
measure x_anc[1] -> x_c[1];
h x_anc[2];
cx x_anc[2],q[1];
cx x_anc[2],a10[0];
cx x_anc[2],q[6];
cx x_anc[2],q[7];
cx x_anc[2],a10[0];
measure a10[0] -> c10[0];
cx x_anc[2],q[15];
h x_anc[2];
measure x_anc[2] -> x_c[2];
h x_anc[3];
cx x_anc[3],q[0];
cx x_anc[3],a11[0];
cx x_anc[3],q[3];
cx x_anc[3],a11[1];
cx x_anc[3],q[10];
cx x_anc[3],q[11];
cx x_anc[3],a11[0];
measure a11[0] -> c11[0];
cx x_anc[3],q[12];
cx x_anc[3],a11[1];
measure a11[1] -> c11[1];
cx x_anc[3],q[14];
h x_anc[3];
measure x_anc[3] -> x_c[3];
