cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENTIONS OFF)

# Define a "project", providing a description and a programming language
project(ex09
    VERSION 1.0
    DESCRIPTION "Exercise 09"
    LANGUAGES CXX
)

find_package(MPI REQUIRED)

message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

# Our "project" will be able to build an executable out of a C++ source file
add_executable(main
    src/main.cpp
)

add_executable(main-parallel
    src/main-parallel.cpp
)

# Our project will be able to build a library
add_library(myLibrary
    src/random.cpp
    src/myStatFunc.cpp
    src/TSP.cpp
)

# This is needed if we keep .h files in the "include" directory
target_include_directories(myLibrary PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    myLibrary PUBLIC 
    MPI::MPI_CXX
)

# Target commands tells which dependencies targer exec have
target_link_libraries(
    main PUBLIC 
    myLibrary
)

target_link_libraries(
    main-parallel PUBLIC 
    myLibrary
    MPI::MPI_CXX
)

# Force the compiler to use the C++17 standard
target_compile_features(main PUBLIC cxx_std_17)
target_compile_features(main-parallel PUBLIC cxx_std_17)
