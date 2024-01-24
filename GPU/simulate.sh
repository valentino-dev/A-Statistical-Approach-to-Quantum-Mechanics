#!/bin/bash

nvcc src/simulation.cu -o "bin/simulation.out"
nvprof bin/simulation.out
