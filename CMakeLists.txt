cmake_minimum_required(VERSION 3.16)
project(IsingModel)

set(CMAKE_CXX_STANDARD_17)

# add_executable(IsingMarkovChainMonteCarloTest ./source/isingMarkovChainMonteCarlo.cpp)
# target_include_directories(IsingMarkovChainMonteCarloTest PRIVATE extern/pybind11/include/)

# pybind11
add_subdirectory(extern/pybind11)
pybind11_add_module(IsingMarkovChainMonteCarlo ./source/isingMarkovChainMonteCarlo.cpp)

add_subdirectory(extern/tracy)
target_link_libraries(IsingMarkovChainMonteCarlo PRIVATE TracyClient)
# target_compile_definitions(IsingMarkovChainMonteCarlo PRIVATE TRACY_ENABLE)

