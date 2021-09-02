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
echo -e "Generating input files for equilibration..."
echo -e "Changing directory to md-configurations..."
cd md-configurations
cp input.dummy input.gas
sed -i 's/temp/1.2/g' input.gas
sed -i 's/npart/108/g' input.gas
sed -i 's/rho/0.05/g' input.gas
sed -i 's/rcut/5.0/g' input.gas
sed -i 's/delta/7.5/g' input.gas
sed -i 's/nstep/10000/g' input.gas
sed -i 's/n_blocks/1/g' input.gas

cp input.dummy input.liquid
sed -i 's/temp/1.1/g' input.liquid
sed -i 's/npart/108/g' input.liquid
sed -i 's/rho/0.8/g' input.liquid
sed -i 's/rcut/2.5/g' input.liquid
sed -i 's/delta/0.2/g' input.liquid
sed -i 's/nstep/10000/g' input.liquid
sed -i 's/n_blocks/1/g' input.liquid

cp input.dummy input.solid
sed -i 's/temp/0.8/g' input.solid
sed -i 's/npart/108/g' input.solid
sed -i 's/rho/1.1/g' input.solid
sed -i 's/rcut/2.2/g' input.solid
sed -i 's/delta/0.12/g' input.solid
sed -i 's/nstep/10000/g' input.solid
sed -i 's/n_blocks/1/g' input.solid

echo -e "Finished generating input files..."
echo -e "Exiting directory..."
cd ..

echo -e "Changing directory to build dir..."
cd build/

echo -e "Executing equilibration..."
for state in solid liquid gas
do
  mkdir -p ex07.1/$state
  mkdir -p ex07.1/$state/results
  mkdir -p ex07.1/$state/frames
  cp ../Primes ex07.1/$state/Primes
  cp ../seed.in ex07.1/$state/seed.in
  cp ../md-configurations/input.$state ex07.1/$state/input.$state
  cp ../md-configurations/config.0 ex07.1/$state/config.0
  cp montecarlomd ex07.1/$state
  cd ex07.1/$state

  ./montecarlomd equilibration -i input.$state --instant yes
  cp -r results results_eq
  cd ..
  cd ..
done

echo -e "Exiting from build dir..."
cd ..
echo -e "Generating input files for restart..."
echo -e "Changing directory to md-configurations..."
cd md-configurations
cp input.dummy input.gas
sed -i 's/temp/1.2/g' input.gas
sed -i 's/npart/108/g' input.gas
sed -i 's/rho/0.05/g' input.gas
sed -i 's/rcut/5.0/g' input.gas
sed -i 's/delta/7.5/g' input.gas
sed -i 's/nstep/50000/g' input.gas
sed -i 's/n_blocks/1/g' input.gas

cp input.dummy input.liquid
sed -i 's/temp/1.1/g' input.liquid
sed -i 's/npart/108/g' input.liquid
sed -i 's/rho/0.8/g' input.liquid
sed -i 's/rcut/2.5/g' input.liquid
sed -i 's/delta/0.2/g' input.liquid
sed -i 's/nstep/50000/g' input.liquid
sed -i 's/n_blocks/1/g' input.liquid

cp input.dummy input.solid
sed -i 's/temp/0.8/g' input.solid
sed -i 's/npart/108/g' input.solid
sed -i 's/rho/1.1/g' input.solid
sed -i 's/rcut/2.2/g' input.solid
sed -i 's/delta/0.12/g' input.solid
sed -i 's/nstep/50000/g' input.solid
sed -i 's/n_blocks/1/g' input.solid

echo -e "Finished generating input files..."
echo -e "Exiting directory..."
cd ..

echo -e "Changing directory to build dir..."
cd build/

echo -e "Executing simulation..."
for state in solid liquid gas
do
  cp ../md-configurations/input.$state ex07.1/$state/input.$state
  cd ex07.1/$state
  ./montecarlomd restart -i input.$state --instant yes
  cd ..
  cd ..
done

cd ..
echo -e "Finished execution."
