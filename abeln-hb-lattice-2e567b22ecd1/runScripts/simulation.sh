#!/bin/sh

# -i              :pdb input file
# -minT,MaxT      :mininum and maximum temperature in reduced units, in between which parallel tempering is performed
# -S              :number of MC steps per round
# -W              :number of rounds, if parallel tempering is on (default), a trial temperature swap is performed after each round
# -noSwap         :turns off parallel tempering, no temperature swaps are performed (use in grand canonic simulation)
# -lT             :linear temperature distribution over the nodes, otherwise linear on 1/temperature  
# -aa             :the input matrix for interaction energies
# -indxWater      :the index of water in the interaction matrix
# -ebmu           :sets the ideal concentration for free peptide chains (no of chains per lattice site), implies simulation in grad canonic ensemble  
# -native         :sets native structure, so that statistics about the number of native contacts (Nint,Next,Ntot) may be computed
# -designT        :sets variability of sequence in structure, implies design program is used instead of simulation

cd ../output

time mpirun -np 2 -output-filename ptemp.out ../source/hblattice -i ../PDBstructures/Tshape1.pdb -minT 0.1 -maxT 0.6 -S 1000 -W 1000 -aa ../matrices/aa_water2.txt -indxWater 20 -lT -hbondStats -CintStats -clusterMoves> ../output/output.out 

cd -
