#!/bin/bash

nvcc src/main.cu
nvprof ./a.out
