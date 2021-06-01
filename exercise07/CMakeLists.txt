cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENTIONS OFF)

# Define a "project", providing a description and a programming language
project(ex07
    VERSION 1.0
    DESCRIPTION "Exercise 07"
    LANGUAGES CXX
)

# Our "project" will be able to build an executable out of a C++ source file
add_executable(montecarlomd
    src/montecarlomd.cpp
)

add_executable(moleculardynamics
    src/moleculardynamics.cpp
)
# Our project will be able to build a library
add_library(myLibrary
    src/random.cpp
    src/mcmd.cpp
    src/MolDyn_NVE.cpp
    src/myStatFunc.cpp
)

# This is needed if we keep .h files in the "include" directory
target_include_directories(myLibrary PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Target commands tells which dependencies targer exec have
target_link_libraries(
    montecarlomd PUBLIC 
    myLibrary
)

target_link_libraries(
    moleculardynamics PUBLIC 
    myLibrary
)
# Force the compiler to use the C++17 standard
target_compile_features(montecarlomd PUBLIC cxx_std_17)
target_compile_features(moleculardynamics PUBLIC cxx_std_17)