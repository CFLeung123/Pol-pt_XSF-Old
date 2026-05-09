# pt_XSF-Old
ifx compiler version, ifort is now deprecated and was discontinued in October 2024


This code computes various correlation functions (gX, lY, g1, l1, gVt, lVt, and 4‑fermion operators) in lattice QCD with Schrödinger Functional (SF) boundary conditions, at tree level and at one loop. The calculations are performed for Wilson and Luscher–Weisz gauge actions, including improvement terms (clover, boundary counterterms).

## Prerequisites

- **Intel oneAPI** (Fortran compiler and MKL). The code uses standard Fortran 90/95 with some OpenMP parallelisation.  
  Load the environment before compiling or running:
  ```bash
  source /opt/intel/oneapi/setvars.sh
  ```
- **GNU make** – for building the executable.
- **Bash** – for the provided runner script.

## Directory Structure

The source files are located in the `source/` directory. After compilation, the executable and module files are copied to `bin/`, and output data goes to `output/`.

```
./
├── compile.sh                 # Compilation script
├── pt_runner.sh               # Execution script (loop over lattice sizes)
├── sizelist                   # List of lattice sizes (one per line)
├── in.dat.src                 # Input template (parameter file)
└── source/                    # Fortran source files (*.f90)
```

## Compilation

Run the compilation script from the top directory:

```bash
./compile.sh
```

What it does:
1. Sources the Intel oneAPI environment.
2. Enters the `source/` directory and runs `make` to build the executable `pt_simple.out` and module files (`*.mod`).
3. Runs `make -t` to update timestamps (optional).
4. Creates `bin/` and `output/` directories.
5. Copies all `*.mod` and `*.out` files from `source/` to `bin/`.

If the compilation fails, check that:
- The Intel Fortran compiler (`ifort`) is in your `PATH`.
- The `Makefile` inside `source/` is correctly configured (it should compile all `.f90` files and link the executable `pt_simple.out`).

## Running the Code

### 1. Prepare the input template

The file `in.dat.src` contains all input parameters with a placeholder `=CFG=` for the lattice size. Example:

```
! Lattice size
L = =CFG=
! Bare quark mass
m0 = 0.0
! Clover coefficient (tree-level)
csw = 1.0
...
```

Modify the parameters according to your calculation (see **Input parameters** below).  
Keep the `=CFG=` placeholder – it will be replaced by the actual lattice size at runtime.

### 2. Create a list of lattice sizes

Write one lattice size per line in a file named `sizelist`. For example:

```
8
10
12
```

### 3. Run the simulation

```bash
./pt_runner.sh
```

This script:
- Sources the Intel environment.
- Sets `OMP_NUM_THREADS=1` (default) and `OMP_STACKSIZE=2G`.
- For each size listed in `sizelist`:
  - Replaces `=CFG=` with the current size in `in.dat.src` and writes `bin/in.dat`.
  - Changes to the `bin/` directory.
  - Executes `pt_simple.out` (timed) and redirects the output to `output.dat`.
  - Moves `output.dat` to `output/<size>.dat`.
  - Returns to the top directory.

After a successful run, the output files will be in `output/` (e.g., `output/8.dat`, `output/10.dat`, …).

## Input Parameters (in.dat.src)

The input file is read by the module `input.f90`. Below are the key parameters (refer to the source for a complete list):

| Parameter         | Description |
|-------------------|-------------|
| `L`               | Lattice size (temporal extent) – replaced by `=CFG=` |
| `m0`              | Bare quark mass (tree level) |
| `csw`             | Clover coefficient (tree level) |
| `ds`              | `ds` boundary counterterm (tree level) |
| `zf`              | `zf` boundary counterterm (tree level) |
| `theta`           | Theta angle for the spatial boundary conditions |
| `G_act`           | Gauge action: `1` = plaquette, `2` = Luscher–Weisz |
| `lambda_gf`       | Gauge fixing parameter |
| `diag_g`          | `1` = take advantage of diagonal gauge propagator, `0` = full |
| `deriv`           | `1` = compute temporal derivative of observable, `0` = no |
| `comp_gX`, `comp_lY`, … | Flags to select which correlation functions to compute (`1` = yes) |
| `comp_tree`, `comp_1loop`, `comp_count` | Compute tree level, 1‑loop diagrams, or counterterms |
| `flav1`–`flav4`   | Flavour indices (1=up, 2=down) |
| `gam_i`           | Gamma matrix index for bilinears (gX, lY) |
| `gam4f_pair`      | Pair of gamma matrices for 4‑fermion operators |
| `kdir`            | Spatial direction for Gk1/Gk2 observables (1,2,3) |

Check the source files (`input.f90`, `parameters.f90`) for the exact meaning of each index and the available gamma matrices.

## Output

For each lattice size the program prints to standard output (which is redirected to `output/<size>.dat`).  
The output contains the values of the requested correlation functions at each time slice. For example, for gX correlators you will see lines like:

```
gAuu_0: 0  0.123456  0.000000
gAuu_1a: 0  -0.001234  0.002345
...
```

Each line typically includes the correlation function name, the time slice index, and the real/imaginary parts of the result.

