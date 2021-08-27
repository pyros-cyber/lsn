#!/bin/bash

echo -e "Checking if build directory exists; if no, run previous equilibration.sh"
if [[ ! -d build ]]
then
	./equilibration.sh
	echo -e "Finish equilibration"
fi

echo -e "Changing directory to build dir..."
cd build/

echo -e "Checking if ex06.1 directory exists; if yes, delete it..."

if [[ -d ex06.1 ]]
then
	rm -rf ex06.1 
	echo -e "ex06.1 directory deleted"
fi

mkdir -p ex06.1

for state in metropolis gibbs
do
  for h in 0 0.02
  do
    mkdir -p ex06.1/$state
    mkdir -p ex06.1/$state/$h
    mkdir -p ex06.1/$state/${h}/results

    cp ../launcher.py ex06.1/launcher.py
    cp ../Primes ex06.1/$state/${h}/Primes
    cp ../seed.in ex06.1/$state/${h}/seed.in
    cp isingmodel ex06.1/$state/${h}/isingmodel

  done
done

echo -e "Changing directory to ex06.1..."
cd ex06.1
echo -e "Launching simulation..."
python3 launcher.py
cd ../../
echo -e "Finished execution"
