cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENTIONS OFF)

# Define a "project", providing a description and a programming language
project(ex03
    VERSION 1.0
    DESCRIPTION "Exercise 03"
    LANGUAGES CXX
)

# set PKG_CONFIG_PATH to include lib directory in my conda env
# in order to be found by pkg-config
#set(ENV{PKG_CONFIG_PATH} "$ENV{CONDA_PREFIX}/lib/pkgconfig/:$ENV{PKG_CONFIG_PATH}")

#enable_testing()

# Our "project" will be able to build an executable out of a C++ source file
add_executable(ex03.1
    src/ex03.1.cpp
)

# Our project will be able to build a library
add_library(myLibrary
    src/random.cpp
    src/myStatFunc.cpp
    src/randomWalk.cpp
)

# This is needed if we keep .h files in the "include" directory
target_include_directories(myLibrary PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Target commands tells which dependencies targer exec have
target_link_libraries(
    ex03.1 PUBLIC 
    myLibrary
)

# Force the compiler to use the C++17 standard
target_compile_features(ex03.1 PUBLIC cxx_std_17)