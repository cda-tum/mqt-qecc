#!/bin/bash
path=$1

gap -q << EOF > $path/info.txt
LoadPackage("QDistRnd");;
hx:=ReadMTXE("$path/hx.mtx");;

hz:=ReadMTXE("$path/hz.mtx");;

Print("dZ = ", DistRandCSS(hx[3], hz[3], 150, 0,2), "\n");
Print("dX = ", DistRandCSS(hz[3], hx[3], 150, 0,2), "\n");

EOF
