cmake_minimum_required(VERSION 3.15)
project(SimpleMolecularDynamics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 REQUIRED)

pybind11_add_module(simple_md
    src/cpp/bindings.cpp
    src/cpp/Particle.cpp
    src/cpp/MDSimulation.cpp
)

target_include_directories(simple_md PRIVATE src/cpp)

set_target_properties(simple_md PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/compiled"
)
