# IsingModel
Solving the Ising model with a markov chain Monte Carlo simulation.

The bulk of the work is stored in the jupyter notebook `montecarlomethod.ipynb`. Here I describe the Ising model and the approach we are going to use to solve it.

I also wrote a version of the Markov Chain Monte Carlo in C++ and used pybind11 
to turn it into a Python module `IsingMarkovChainMonteCarlo.so`.

The file `meanfieldapprox.py` is a script I wrote to plot the solution to the mean field approximation of the Ising Model.