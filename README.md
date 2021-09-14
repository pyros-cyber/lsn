# lsn
This repository contains the exercises for the Lab of Numeric Simulation, A.A. 2020/21.

## Installation ##
To download the exercises the easiest way is cloning the repository:
```
$ git clone https://github.com/pyros-cyber/lsn.git
```

### Dependencies ###
In order to compile and run the C++ programs and to visualize the results in the Jupyter-Notebooks there are many dependencies; among the most importants:
- C++ Compiler (at least C++17)
- [CMake](https://cmake.org/) (at least version 3.16)
- [Python](https://www.python.org/downloads/) (at least version 3.7.x)
- [Root](https://root.cern/) (latest version)
- [Open-MPI](https://www.open-mpi.org/) (latest version)

If you're willing to use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html) (or if you already do), all the dependencies can be installed via:
```
$ conda env create -n <envname> -f requirements.yml
```

### Directories structure ###
The general structure of each C++ project directory is, e.g for `exercise01`:
```
$ tree exercise01/
exercise01
├── CMakeLists.txt
├── exercise01.ipynb
├── include
│   ├── myStatFunc.h
│   └── random.h
├── LSN_Exercises_01.ipynb
├── Primes
├── seed.in
└── src
    ├── ex01.1.cpp
    ├── ex01.2.cpp
    ├── ex01.3.cpp
    ├── myStatFunc.cpp
    └── random.cpp
```

## Building and executing ##
If all dependecies have been met, to generate the results (this step can be skipped, as in the jupyter-notebooks there are already the results obtained from previous runs, but if you want to make sure they are correct keep on reading), for instance for `exercise01`, do:
```
$ cd exercise01/
($ conda activate <envname>)
$ cmake -S . -B build/ -DCMAKE_BUILD_TYPE=Release
$ cmake --build build
$ cd build
$ ./<executables-files>
```
For certain exercises bash scripts that compile and run the program are provided, in order to simplify the input paramaters generations and feed them into the executable; more detailed instructions per exercise can be found in the related notebooks (`exerciseNN.ipynb`).

> The CLI interface of all the exercises has been built using [Taywee/args](https://github.com/Taywee/args)

## License ##
This repository is licensed under the terms of the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html) and is available for free. For more details see [here](https://github.com/pyros-cyber/lsn/blob/main/LICENSE.md).

