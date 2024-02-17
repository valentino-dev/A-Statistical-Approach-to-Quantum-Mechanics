#!/bin/bash

nvcc src/simulation.cu -o "bin/simulation.out"
a=0.5
m0=0.5
lambda=0.0
mu=2.0
runs=2
n=15
SII=10
RS=50
TpB=512
MCS=136
MCI=100
xjpT=2
nvprof bin/simulation.out "data/data_fig6.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "0.0"
