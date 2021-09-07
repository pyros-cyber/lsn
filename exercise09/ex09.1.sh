#!/bin/bash

echo -e "Checking if build directory exists; if yes, delete it..."

if [[ -d build ]]
then
	rm -rf build
	echo -e "Build directory deleted."
fi

echo -e "Compiling program.."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
echo -e "and building.."
cmake --build build

echo -e "Changing directory to build dir..."
cd build/

cp ../Primes .
cp ../seed.in .

if [[ -d results ]]
then
	rm -rf results
fi

mkdir -p results

./main --Npop 1000 --Ngen 150 --shape square
./main --Npop 1000 --Ngen 150 --shape circumference
