# A Statistical Approach to Quantum Mechanics

Monte Carlo simulations and analysis for a statistical approach to quantum mechanics. The repository contains a CUDA implementation, an earlier CPU implementation, analysis scripts, and existing data and plots. The code is exploratory research software; a successful run alone does not establish the validity of its results.

## Repository at a glance

| Path | Contents |
| --- | --- |
| `src/simulation/simulation.cu` | CUDA simulation |
| `CPU/Simulation/Simulation/Simulation.cpp` | Original CPU simulation |
| `src/processing/` | Python analysis scripts |
| `bash/` | Example simulation and profiling scripts |
| `data/`, `plots/` | Existing datasets and figures |
| `literature/` | Reference material |
| `src/legacy/` | Earlier simulation and analysis code |

## CUDA quick start

You need an NVIDIA GPU and the CUDA toolkit, including `nvcc` and cuRAND. From the repository root:

```bash
mkdir -p bin data
nvcc src/simulation/simulation.cu -o bin/simulation.out
bin/simulation.out data/example.csv 1 1 2 2 2 32 1 4 0.5 0.5 0.0 2.0 0.0
```

This small smoke run writes four configurations, each with 64 semicolon-separated site values and no header. It is intended to check execution, not to produce equilibrated scientific results. The output file is overwritten if it already exists.

The program expects all 14 arguments, in this order:

```text
output_path runs SII ntherm xjPerThread nhit threadsPerBlock MCS MCI a m0 lambda mu_sq f_sq
```

`SII` is the number of iterations between measurements, `ntherm` is the initial thermalization count, `nhit` is the number of update attempts per site, and `MCS` and `MCI` control parallel configurations and iterations. The number of sites per configuration is `xjPerThread * threadsPerBlock`. See the source and scripts for the remaining simulation parameters. Always provide the full argument list: the fallback settings do not safely supply an output path. Although `f_sq` is accepted, the active simulation uses `potential_1`, which does not use it.

The scripts in `bash/` demonstrate larger parameter sets and should be run from the repository root. Most invoke `nvprof`, which must be installed separately and may not support newer CUDA environments. `bash/test_simulation.sh` is a profiling run, not an automated correctness test. `bash/process.sh` refers to a CUDA processing file that is no longer present.

## CPU simulation

The CPU program's settings are defined in its source. It can be built with a C++ compiler:

```bash
g++ -std=c++17 CPU/Simulation/Simulation/Simulation.cpp -o /tmp/quantum-cpu
```

It writes `data_fig6.csv` in its working directory. **The current implementation has an out-of-bounds read at the first site**, detected by AddressSanitizer; do not rely on its results until the boundary handling is corrected and validated.

## Analysis

The Python scripts use NumPy, pandas, SciPy, Matplotlib and Dask; plotting is configured to use LaTeX. Install these dependencies in a Python environment and review the input paths and selected functions in the scripts before running them. For example, `processing.py` currently reads `data/correlation.csv` by default. `processingPaper.py` also imports `olib` from the legacy directory and contains machine-specific paths, so it needs local configuration before use.

## Verification status

- A small CUDA run compiled and completed on an NVIDIA GPU, producing four rows of 64 finite values. This is a smoke test, not a statistical or numerical validation. The program may print an infinite time estimate at the start of a short run.
- The CPU program compiles, but AddressSanitizer reports a heap-buffer-overflow in its action calculation.
- Python files compile syntactically. The analysis was not run in the checked environment because its scientific Python dependencies were incomplete.

Existing plots and datasets are included for reference; reproduce and independently validate results before drawing conclusions.
