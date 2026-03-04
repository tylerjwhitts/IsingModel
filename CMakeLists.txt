cmake_minimum_required(VERSION 3.16)
project(IsingModel)

set(CMAKE_CXX_STANDARD_17)

# pybind11
add_subdirectory(./pybind11)
pybind11_add_module(IsingMarkovChainMonteCarlo ./source/isingMarkovChainMonteCarlo.cpp)