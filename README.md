# IsingModel
Solving the Ising model with a markov chain Monte Carlo simulation.

The bulk of the work is stored in the jupyter notebook `montecarlomethod.ipynb`. Here I describe the Ising model and the approach we are going to use to solve it.

I also wrote a version of the Markov Chain Monte Carlo in C++ and used pybind11 
to turn it into a Python module `IsingMarkovChainMonteCarlo.so`.

The file `meanfieldapprox.py` is a script I wrote to plot the solution to the mean field approximation of the Ising Model.

## Compiling the C++ Markov Chain Library
1. Create a `build` directory. 
2. Enter the directory `cd build`.
3. Generate the build files using `cmake ..`.
4. Build the project using `make`. This will produce a file called `IsingMarkovChainMonteCarlo.cpython-314-darwin.so` in the `build` directory.
5. To make this library available to the Jupyter notebook move it to the main project directory:
```shell
mv IsingMarkovChainMonteCarlo.cpython-314-darwin.so ../IsingMarkovChainMonteCarlo.so
```
6. You can now import the library as you would any other Python module.

NOTE: The CMake process also creates a test executable for the C++ code. In order for this to compile correctly the following header directories must be in your include compilation flags:
- The tracy directory: `${workspaceFolder}/extern/tracy/public/`
- The pybind11 directory: `${workspaceFolder}/extern/pybind11/include/`
- The directory where your `Python.h` file is stored. For a homebrew installed python this will look something like `/opt/homebrew/Cellar/python@3.14/3.14.2_1/Frameworks/Python.framework/Versions/3.14/include/python3.14/`.