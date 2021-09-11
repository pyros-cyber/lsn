#!/bin/bash

echo -e "Changing directory to build dir..."
cd build/

mkdir ex10.2

cp ../Primes ex10.2/
cp ../seed.in ex10.2/
cp main-parallel ex10.2/

cd ex10.2

if [[ -d results ]]
then
	rm -rf results
fi

mkdir -p results

mpiexec -n 4 main-parallel --Npop 1000 --Ngen 150 --shape square
mpiexec -n 4 main-parallel --Npop 1000 --Ngen 150 --shape circumference

cd ../..
echo "Finish execution"
