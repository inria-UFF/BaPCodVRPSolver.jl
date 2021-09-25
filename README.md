# BaPCodVRPSolver

The BaPCodVRPSolver.jl package is a Julia interface to use the [VRPSolver](https://vrpsolver.math.u-bordeaux.fr/) in Linux operating systems. 

## Requirements

- Linux (preferably Ubuntu)
- [Julia](https://julialang.org/downloads/oldreleases/) between 1.2 and 1.5 (preferably 1.4)
- [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) 12.9 or 12.10 

Other versions of Julia and CPLEX are not guaranteed to work.

## Installation

Open the Julia REPL and type:
```
    ]add https://github.com/inria-UFF/BaPCodVRPSolver.jl.git
```

Set the *LD_LIBRARY_PATH* environment variable with the absolute path to the subdirectory of your CPLEX installation which contains the executable files and shared libraries.
For example, if your CPLEX is installed at */opt/ibm/ILOG/CPLEX_Studio1210* and you are using Bash, you can declare it in the ~/.bashrc:

```
export LD_LIBRARY_PATH=/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/:$LD_LIBRARY_PATH
```

Next, set the *BAPCOD_RCSP_LIB* environment variable with the absolute path to the *BaPCod_RCSP.so* library (requested through the form on the [home page](https://vrpsolver.math.u-bordeaux.fr/)).
For example, if the *BaPCod_RCSP.so* is at */path/to/lib/BaPCod_RCSP.so* and you are using Bash, you can declare it in the ~/.bashrc:

```
export BAPCOD_RCSP_LIB=/path/to/lib/BaPCod_RCSP.so
```

## Running a application

All the demos available in the [VRPSolver website](https://vrpsolver.math.u-bordeaux.fr/) will work after replacing `using VrpSolver` with `using BaPCodVRPSolver` in the file *src/run.jl*.

For example, the [CVRP demo](https://vrpsolver.math.u-bordeaux.fr/cvrpdemo.zip) can be invoked for the instance X-n101-k25 using an upper bound of 27591.1 as follows:

```
julia src/run.jl data/X/X-n101-k25.vrp -u 27591.1
```

