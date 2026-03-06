cmake_minimum_required(VERSION 3.16)
project(IsingModel)

set(CMAKE_CXX_STANDARD 17)

# pybind11 from submodule
add_subdirectory(extern/pybind11)

# tracy
add_subdirectory(extern/tracy)

# Create a python module
pybind11_add_module(IsingMarkovChainMonteCarlo ./source/isingMarkovChainMonteCarlo.cpp)

target_link_libraries(IsingMarkovChainMonteCarlo PRIVATE Tracy::TracyClient)
target_compile_definitions(IsingMarkovChainMonteCarlo PRIVATE TRACY_ENABLE)

# Testing the original C++ program
add_executable(IsingMarkovChainMonteCarloTest ./source/isingMarkovChainMonteCarlo.cpp)
target_link_libraries(IsingMarkovChainMonteCarloTest PRIVATE pybind11::embed Tracy::TracyClient)
target_compile_definitions(IsingMarkovChainMonteCarloTest PRIVATE TRACY_ENABLE)
