#!/bin/bash
# Bash script to compile and execute program

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
echo -e "Generating input files..."
echo -e "Changing directory to md-configurations..."
cd md-configurations
cp input.dummy input.gas
sed -i 's/temp/1.2/g' input.gas
sed -i 's/npart/108/g' input.gas
sed -i 's/rho/0.05/g' input.gas
sed -i 's/rcut/5.0/g' input.gas
sed -i 's/delta/0.0005/g' input.gas
sed -i 's/nstep/1000/g' input.gas
sed -i 's/iprint/1000/g' input.gas
sed -i 's/measure_time_interval/10/g' input.gas
sed -i 's/n_blocks/1/g' input.gas

cp input.dummy input.liquid
sed -i 's/temp/1.1/g' input.liquid
sed -i 's/npart/108/g' input.liquid
sed -i 's/rho/0.8/g' input.liquid
sed -i 's/rcut/2.5/g' input.liquid
sed -i 's/delta/0.0005/g' input.liquid
sed -i 's/nstep/1000/g' input.liquid
sed -i 's/iprint/1000/g' input.liquid
sed -i 's/measure_time_interval/10/g' input.liquid
sed -i 's/n_blocks/1/g' input.liquid

cp input.dummy input.solid
sed -i 's/temp/0.8/g' input.solid
sed -i 's/npart/108/g' input.solid
sed -i 's/rho/1.1/g' input.solid
sed -i 's/rcut/2.2/g' input.solid
sed -i 's/delta/0.0005/g' input.solid
sed -i 's/nstep/1000/g' input.solid
sed -i 's/iprint/1000/g' input.solid
sed -i 's/measure_time_interval/10/g' input.solid
sed -i 's/n_blocks/1/g' input.solid

echo -e "Finished generating input files..."
echo -e "Exiting directory..."
cd ..

echo -e "Changing directory to build dir..."
cd build/

echo -e "Executing program..."
for state in solid liquid gas
do
  mkdir -p ex04.1/$state
  mkdir -p ex04.1/$state/results
  mkdir -p ex04.1/$state/frames
  cp ../Primes ex04.1/$state/Primes
  cp ../seed.in ex04.1/$state/seed.in
  cp ../md-configurations/input.$state ex04.1/$state/input.$state
  cp ../md-configurations/config.0 ex04.1/$state/config.0
  cp moleculardynamics ex04.1/$state
  cp ../md-configurations/clean.sh ex04.1/$state
  cd ex04.1/$state
  for ((irnd=1; irnd<7; ++irnd))
  do
    if [[ "$irnd" -eq "1" ]]
    then
      ./moleculardynamics equilibration -i input.$state
    else
      ./moleculardynamics  restart -i input.$state
    fi
    cp -r results results$irnd
  done
  cd ..
  cd ..
done
echo -e "Finished execution."
