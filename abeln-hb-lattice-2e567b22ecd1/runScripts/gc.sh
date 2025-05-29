#!/bin/sh

cd ../output

time mpirun -np 1 -output-filename ptemp.out ../source/hblattice -i ../PDBstructures/T3LT.pdb -minT 0.10 -maxT 0.10 -S 1000 -W 1000 -aa ../matrices/aa_water2.txt -indxWater 20 -lT -noSwap -ebmu 0.003 > output.out 


cd -

