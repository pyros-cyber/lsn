#!/bin/bash

echo -e "Changing directory to build dir..."
cd build/

echo "Generating parameters..."

mkdir ex08.3

cp qmc1d ex08.3/
cp qmc1d_const ex08.3/
cp -r ../input/ ex08.3/
cp ../input.dat ex08.3/

cd ex08.3
mkdir results

#PIGS
echo "Executing PIGS algorithm with variational trial wave function"
for tau in 05 1 2 3
do
	cp input/input_t$tau.pigs input.dat
	./qmc1d 	#variational trial wave function

	for res in probability potential kinetic
	do
		cp $res.dat results/$res.var_t$tau.dat
	done;echo;echo

done


echo "Executing PIGS algorithm with constant trial wave function"
for tau in 3 5 6
do
	cp input/input_t$tau.pigs input.dat
	./qmc1d_const 	#constant trial wave function

	for res in probability potential kinetic
	do
		cp $res.dat results/$res.const_t$tau.dat
	done;echo;echo

done


#PIMC
echo "Executing PIMC algorithm with variational trial wave function"
for temp in 0_25 1_25 5 50
do
	cp input/input_$temp.pimc input.dat
	./qmc1d

	for res in probability potential kinetic
	do
		cp $res.dat results/$res.T$temp.dat
	done;echo;echo

done

for res in probability potential kinetic
do
	rm -rf $res.dat
done

cd ../..
echo "Finished execution."
