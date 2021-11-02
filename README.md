# BaPCodVRPSolver

The BaPCodVRPSolver.jl package is a Julia interface to use the [VRPSolver](https://vrpsolver.math.u-bordeaux.fr/) in
Linux operating systems.

This package is *only for academic use*.

## Requirements

- [Julia](https://julialang.org/downloads/oldreleases/) versions 1.0 -- 1.5.4 (versions 1.6 and later will be supported
  in the future)
- [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) versions 12.9 and higher

## Installation

Open the Julia REPL and type:
```
    ]add https://github.com/inria-UFF/BaPCodVRPSolver.jl.git
```

On Linux, set the *LD_LIBRARY_PATH* environment variable with the absolute path to the subdirectory of your CPLEX
installation which contains the executable files and shared libraries.  For example, if your CPLEX is installed at
*/opt/ibm/ILOG/CPLEX_Studio1210* and you are using Bash, you can declare it in the ~/.bashrc:

```
export LD_LIBRARY_PATH=/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/:$LD_LIBRARY_PATH
```

On Windows, be sure that the *PATH* environment variable contains the folder with CPLEX dynamic library.

Next, set the *BAPCOD_RCSP_LIB* environment variable with the absolute path to the BaPCod shared library (which has
*.so* extension on linux, *.dylib* extension on Mac OS, and *.dll* extension on Windows).
For example, if you are using Bash on Linux, you can declare it in the ~/.bashrc:

```
export BAPCOD_RCSP_LIB=/path/to/lib/libbapcod-shared.so
```

## Producing BaPCod shared library

If the BaPCod shared library you have does not work for you or you do not have one, you can produce it in the following
way. Download BaPCod source code on its web-page: https://bapcod.math.u-bordeaux.fr/. Then request the BCP_RCSP compiled
library from Ruslan(point)Sadykov(at)inria(point)fr. Install BaPCod together with the BCP_RCSP library it using installations
instruction in README.md file or in the user guide (available on its web-page).

Then run the following command from BapcodFramework folder 

```
cmake --build build --config Release --target bapcod-shared
```

This will produce shared library file build/Bapcod/libbapcod-shared.so on Linux, build/Bapcod/libbapcod-shared.dylib on Mac OS,
and build/Bapcod/Release/bapcod-shared.dll on Windows.


## Running a application

Firstly, add the folowing dependences:

```
   ]add JuMP, CPLEX, ArgParse
```

All the demos available in the [VRPSolver website](https://vrpsolver.math.u-bordeaux.fr/) will work after replacing `using VrpSolver` with `using BaPCodVRPSolver` in the file *src/run.jl*.

For example, the [CVRP demo](https://vrpsolver.math.u-bordeaux.fr/cvrpdemo.zip) can be invoked (after making the aforementioned replacement) for the instance X-n101-k25 using an upper bound of 27591.1 as follows:

```
julia src/run.jl data/X/X-n101-k25.vrp -u 27591.1
```

