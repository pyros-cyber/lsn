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

mkdir ex10.1

cp ../Primes ex10.1/
cp ../seed.in ex10.1/
cp main ex10.1/

cd ex10.1

if [[ -d results ]]
then
	rm -rf results
fi

mkdir -p results

./main --Npop 1000 --Ngen 200 --ntemp 500 --nsteps 2000 --shape square
./main --Npop 1000 --Ngen 200 --ntemp 500 --nsteps 100 --shape circumference

cd ../..
echo "Finish execution."
