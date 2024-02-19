#!/bin/bash

nvcc src/processing.cu -o "bin/processing.out"
nvprof bin/processing.out
