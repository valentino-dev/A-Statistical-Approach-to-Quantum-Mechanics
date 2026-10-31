# A Statistical Approach to Quantum Mechanics

The project contains CPU and GPU implementations in the same branch:

- `CPU/Simulation/Simulation/Simulation.cpp`: original CPU simulation. Its settings are defined in the source; it writes `data_fig6.csv` in the working directory.
- `src/simulation/simulation.cu`: current CUDA simulation. The scripts in `bash/` contain example parameters and should be run from the repository root. They require the CUDA toolkit; most also invoke `nvprof`.
- `src/processing/`: GPU branch's Python analysis scripts. These use scientific Python packages, including NumPy, pandas, SciPy, Matplotlib and Dask; plots may require a LaTeX installation.
- `src/legacy/`: earlier simulation and analysis files retained from main. They are historical and have not been reconciled with the newer analysis.
- `data/`, `plots/`, and `literatur/`: existing data, results and references from the project.

For example, compile the current GPU simulation from the repository root with `nvcc src/simulation/simulation.cu -o bin/simulation.out` (create `bin/` first if needed). See the simulation scripts for arguments. `bash/process.sh` is an old script referring to a CUDA processing file that is no longer present; use the Python analysis scripts instead after reviewing their configured input paths.

Make new changes on short-lived branches based on `main`, then merge them back so CPU, GPU and analysis changes stay together. The older hardware branches remain available for their history.
