#!/bin/bash

nvcc src/simulation.cu -o "bin/simulation.out"
a=0.5
m0=0.5
lambda=0.0
mu=2.0
runs=1
n=10
SII=1
RS=0
TpB=512
MCS=1
MCI=10000
xjpT=2
f=0.0
nvprof bin/simulation.out "data/data_fig5.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "${f}"
