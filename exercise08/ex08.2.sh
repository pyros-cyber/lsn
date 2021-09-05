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

echo "Generating parameters..."

mkdir ex08.2

cp ../Primes ex08.2/Primes
cp ../seed.in ex08.2/seed.in
cp main ex08.2/

cd ex08.2

if [ -f variational.dat ]; then
    rm variational.dat
fi

touch variational.dat

MU_MIN=0.7
MU_MAX=0.9
MU_STEP=0.01

SIGMA_MIN=0.5
SIGMA_MAX=0.8
SIGMA_STEP=0.01

for mu in `seq $MU_MIN $MU_STEP $MU_MAX`
do
  for sigma in `seq $SIGMA_MIN $SIGMA_STEP $SIGMA_MAX`
  do
    echo "$mu $sigma"
    echo ${mu/,/.}
    echo ${sigma/,/.}
    ./main  -m ${mu/,/.} -s ${sigma/,/.}
    tail -1 energy.dat | awk -v mu="${mu/,/.}" -v sigma="${sigma/,/.}" '{print mu, sigma, $2, $3}' >> variational.dat
  done
done

RESULTS=($(awk '{print $3, $1, $2}' variational.dat | sort --numeric-sort | head -1))

echo "Optimal parameters are: mu=${RESULTS[1]}, sigma=${RESULTS[2]}"

echo "Results saved in 'variational.dat'..."
echo "All done!"

echo "Now we run a simulation with optimal parameters and 10 000 000 steps"

./main -m ${RESULTS[1]} -s ${RESULTS[2]} -M 10000000

cd ../..
echo "Finished execution."
