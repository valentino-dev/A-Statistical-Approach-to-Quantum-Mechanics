#!/bin/bash

nvcc src/simulation.cu -o "bin/simulation.out"

bin/simulation.out "data/data_fig7a.csv" "2" "10" "50" "2" "15" "512" "136" "100" "0.25" "0.5" "1.0" "2.0" "0.0"
bin/simulation.out "data/data_fig7b.csv" "2" "10" "50" "2" "15" "512" "136" "100" "0.25" "0.5" "1.0" "2.0" "0.5"
bin/simulation.out "data/data_fig7c.csv" "2" "10" "50" "2" "15" "512" "136" "100" "0.25" "0.5" "1.0" "2.0" "1.0"
bin/simulation.out "data/data_fig7d.csv" "2" "10" "50" "2" "15" "512" "136" "100" "0.25" "0.5" "1.0" "2.0" "1.5"
bin/simulation.out "data/data_fig7e.csv" "2" "10" "50" "2" "15" "512" "136" "100" "0.25" "0.5" "1.0" "2.0" "2.0"
