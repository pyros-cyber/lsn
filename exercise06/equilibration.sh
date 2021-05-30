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

echo -e "Checking if equilibration directory exists; if yes, delete it..."

if [[ -d equilibration ]]
then
	rm -rf equilibration
	echo -e "Equlibration directory deleted"
fi

mkdir -p equilibration

for state in metropolis gibbs
do
  for h in 0 0.02
  do
    for t in 2.0 0.5 0.95
    do
      mkdir -p equilibration/$state
      mkdir -p equilibration/$state/$h
      mkdir -p equilibration/$state/${h}/${t}
      mkdir -p equilibration/$state/${h}/${t}/results

      cp ../Primes equilibration/$state/${h}/${t}/Primes
      cp ../seed.in equilibration/$state/${h}/${t}/seed.in
      cp isingmodel equilibration/$state/${h}/${t}/isingmodel
      cp ../input.dummy equilibration/$state/${h}/${t}/input.dat

      sed -i "s/h_dummy/${h}/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/nblk_dummy/1/g" equilibration/$state/${h}/${t}/input.dat
      sed -i "s/nstep_dummy/5000/g" equilibration/$state/${h}/${t}/input.dat

      #if [ "$state" == "metropolis" ];
      #then
      #  sed -i "s/metro_dummy/1/g" equilibration/$state/${h}/${t}/input.dat
      #elif [ "$state" == "gibbs" ];
      #then
      #  sed -i "s/metro_dummy/0/g" equilibration/$state/${h}/${t}/input.dat
      #fi

      cd equilibration/$state/${h}/${t}
      sed -i "s/temp_dummy/${t}/g" input.dat
      ./isingmodel equilibration -a $state 
      mv results results_eq
      cd ../../../../
      done
  done
done
cd ../
echo -e "Finished execution"
