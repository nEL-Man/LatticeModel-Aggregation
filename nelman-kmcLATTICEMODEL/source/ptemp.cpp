import math
import random
import sys
from mpi4py import MPI
from typing import List, Tuple, Optional

# Constants
MAX_PROCS = 32
MASTER = 0
SWAP_REQUEST = 1
ALL_DONE = 2
LAMBDA_COR = False

# Global variables
betas: List[float] = [0.0] * MAX_PROCS
betaID2node: List[int] = [0] * MAX_PROCS
node2betaID: List[int] = [0] * MAX_PROCS
energies: List[int] = [0] * MAX_PROCS
nChains: List[int] = [0] * MAX_PROCS
swapacc: List[int] = [0] * MAX_PROCS
swaptry: List[int] = [0] * MAX_PROCS

myBeta: float = 0.0
swapT: bool = True

# Assuming these classes are defined elsewhere
class Design:
    pass

class MonteCarlo:
    pass

class App:
    pass

def init_betas(min_beta: float, max_beta: float, num_procs: int, linear_T: bool) -> None:
    if num_procs == 1:
        betas[0] = (min_beta + max_beta) / 2
        betaID2node[0] = 0
        node2betaID[0] = 0
    else:
        if not linear_T:
            step = (max_beta - min_beta) / (num_procs - 1)
            for i in range(num_procs):
                betas[i] = min_beta + i * step
                betaID2node[i] = i
                node2betaID[i] = i
                print(f"node {i} beta {betas[i]} T {0.01 / betas[i]}")
        else:
            max_T = 0.01 / min_beta
            min_T = 0.01 / max_beta
            temperatures = [0.0] * num_procs
            step = (max_T - min_T) / (num_procs - 1)
            for i in range(num_procs):
                temperatures[num_procs - i - 1] = min_T + i * step
            
            for i in range(num_procs):
                betas[i] = 0.01 / temperatures[i]
                betaID2node[i] = i
                node2betaID[i] = i
                print(f"node {i} beta {betas[i]} T {temperatures[i]}")

def try_swap(betaID1: int, betaID2: int) -> bool:
    print(f"trial: {betaID1} {betaID2}")
    node1 = betaID2node[betaID1]
    node2 = betaID2node[betaID2]
    beta1 = betas[betaID1]
    beta2 = betas[betaID2]
    real_beta1 = 100 * beta1
    real_beta2 = 100 * beta2

    energy1 = energies[node1]
    energy2 = energies[node2]

    nChains1 = nChains[node1]
    nChains2 = nChains[node2]

    dE = float(energy1 - energy2)
    dB = beta1 - beta2
    
    if LAMBDA_COR:
        prefactor = pow(real_beta1 / real_beta2, (3/2) * (nChains1 - nChains2))
    else:
        prefactor = 1.0

    p_acc_swap = prefactor * math.exp(dE * dB)
    swaptry[betaID1] += 1
    swaptry[betaID2] += 1

    if random.random() < p_acc_swap:
        print(f"SWAPPING {beta1} with energy: {energy1} at node {node1} "
              f"with {beta2} with energy: {energy2} at node {node2}")
        betaID2node[betaID1] = node2
        betaID2node[betaID2] = node1
        node2betaID[node1] = betaID2
        node2betaID[node2] = betaID1
        swapacc[betaID1] += 1
        swapacc[betaID2] += 1
        return True
    else:
        return False

def main():
    # Initialize MPI
    comm = MPI.COMM_WORLD
    num_procs = comm.Get_size()
    my_rank = comm.Get_rank()
    name = MPI.Get_processor_name()
    
    App.myRank = my_rank  # Assuming App is a class with this attribute

    # Check number of processes
    if num_procs > MAX_PROCS:
        print(f"Too many processes: {num_procs}")
        sys.exit(1)

    print(f"Process {my_rank} out of {num_procs} realp {name}")

    # Initialize random number generator with rank
    random.seed(my_rank)

    # Initialize App (assuming this is done elsewhere)
    # App.init(argc, argv)  # In Python, we'd handle command line args differently

    # Initialize Monte Carlo and/or Design
    design = None
    mc_sim = None

    if App.designProgram:
        design = Design()
        if not App.useOldFormat:
            design.init(App.fn_in)
        else:
            design.initOldFormat(App.fn_in)
    else:
        mc_sim = App.getMonteCarlo()

    print("Finished initializing")

    # Initialize betas (only on master)
    if my_rank == MASTER:
        init_betas(App.minBeta, App.maxBeta, num_procs, App.linear_T)

    # Broadcast beta information
    betaID2node = comm.bcast(betaID2node if my_rank == MASTER else None, root=MASTER)
    betas = comm.bcast(betas if my_rank == MASTER else None, root=MASTER)
    node2betaID = comm.bcast(node2betaID if my_rank == MASTER else None, root=MASTER)
    
    myBeta = betas[node2betaID[my_rank]]
    print("Betas set")

    # Main simulation loop
    for swaps in range(mc_sim.total_swaps):
        my_beta_changed = True
        my_energy = 0

        if App.designProgram:
            design.designProcedure(myBeta, mc_sim.steps)
        else:
            my_energy = mc_sim.MC(myBeta, node2betaID[my_rank], my_beta_changed)

        # Gather energies and chain counts
        energies = comm.gather(my_energy, root=MASTER)
        my_n_chains = mc_sim.grancan.getNumFreeChains()
        nChains = comm.gather(my_n_chains, root=MASTER)

        if swapT and my_rank == MASTER:
            start = 0 if random.random() < 0.5 else 1  # swap even or uneven
            for i in range(start, num_procs - 1, 2):
                try_swap(i, i + 1)

        # Broadcast updated beta information
        node2betaID = comm.bcast(node2betaID if my_rank == MASTER else None, root=MASTER)
        betaID2node = comm.bcast(betaID2node if my_rank == MASTER else None, root=MASTER)
        myBeta = betas[node2betaID[my_rank]]

    # Finalize
    if App.designProgram:
        design.finalize()
    else:
        mc_sim.getMCStats(num_procs, my_rank)

    if my_rank == MASTER:
        for i in range(num_procs):
            ratio = swapacc[i] / swaptry[i] if swaptry[i] > 0 else 0.0
            print(f"beta {i} trials={swaptry[i]} acc={swapacc[i]} ratio={ratio}")

    print(f"Finalizing rank {my_rank}")
    MPI.Finalize()
    print(f"Close down process {my_rank}")
    sys.exit(0)

if __name__ == "__main__":
    main()
