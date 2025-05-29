#!/bin/sh

cd ../output

time mpirun -np 1 ../source/hblattice -i ../PDBstructures/Tshape.pdb -minT 0.1 -maxT 0.1 -S 1000000 -aa ../matrices/aa_water2.txt -indxWater 20 -designT 30  > output.out 

cd -
