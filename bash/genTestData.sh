
#!/bin/bash

nvcc src/simulation.cu -o "bin/simulation.out"
runs=1
SII=1
RS=1000
xjpT=2
n=10
TpB=512
MCS=1
MCI=10000
a=0.5
m0=0.5
lambda=0.0
mu=2.0
f=0.0
nvprof bin/simulation.out "data/data_fig5_10k_RS_1k.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "${f}"
