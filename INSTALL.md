## cmake

Something like the following should work, assuming that you have GMP and GLPK locatable.

```
cmake -B build -S .
cmake --build build
cmake --install --prefix=install
```

You should get an `install` subdirectory with a library and a `scallop` executable under `bin`.

One way to satisfy the dependencies is to install conda/mamba, say from https://github.com/conda-forge/miniforge, and then
```
mamba create -n scallop-build cmake gmp glpk
```

To run `scallop`, you may find that you need to `export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_PREFIX}/lib` so that the executable can find the glpk shared libary.
On MacOS, the environment variable is `DYLD_LIBRARY_PATH`.

### Using gurobi with cmake

Maybe the following will work?  Edit the configure line:
```
cmake -B build -S . -DUSE_GUROBI
```


## With makefiles

To install, you should be able to simply type:

make

Note that scallop requires the development libraries for GLPK and GMP.  
On linux or OSX, these should be available through any package manager 
(say, fink, on OSX).  

Gurobi (http://www.gurobi.com/) is a state-of-the-art proprietary 
linear/integer programming solver.  It is often far faster than GLPK or 
EXLP, and they offer a free academic license.  Scallop supports using 
Gurobi as a solver: to use Gurobi, edit the GURDIR path in the makefile 
to point to your Gurobi directory.  You may also need to edit the 
GURLIB line to make it the appropriate version.  Then to compile, run 

make scallop_with_gurobi

To test, you can run ./scallop -mGUROBI abAB, which should give you 0.5.
