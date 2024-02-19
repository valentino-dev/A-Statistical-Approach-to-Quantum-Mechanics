#!/bin/bash

nvcc src/simulation.cu -o "bin/simulation.out"
a=0.25
m0=0.5
lambda=1.0
mu=2.0
runs=4
n=15
SII=10
RS=50
TpB=512
MCS=136
MCI=100
xjpT=2
nvprof bin/simulation.out "data/data_fig10_0.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "0.0"
nvprof bin/simulation.out "data/data_fig10_1.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "0.5"
nvprof bin/simulation.out "data/data_fig10_2.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "1.0"
nvprof bin/simulation.out "data/data_fig10_3.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "1.5"
nvprof bin/simulation.out "data/data_fig10_4.csv" "${runs}" "${SII}" "${RS}" "${xjpT}" "${n}" "${TpB}" "${MCS}" "${MCI}" "${a}" "${m0}" "${lambda}" "${mu}" "2.0"
