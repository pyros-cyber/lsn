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
echo -e "Changing directory to build dir..."
cd build/
echo -e "Executing program..."
echo -e "Using uniform transition probability"
for cmd in origin far pos
do 
	./h-atom $cmd -s ground -t uniform --step 1.22
	./h-atom $cmd -s first-excited -t uniform --step 2.94
done
echo -e "Using gaussian transition probability"
for cmd in origin far pos
do 
	./h-atom $cmd -s ground -t gaussian --step 1.5
	./h-atom $cmd -s first-excited -t gaussian --step 3.7
done
cd ..
echo -e "Finished execution"
