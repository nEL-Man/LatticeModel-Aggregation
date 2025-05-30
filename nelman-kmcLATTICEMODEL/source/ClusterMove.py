import random
from typing import List, Tuple
from Cluster import Cluster
from Lattice import Lattice
from Chain import Chain
from Residue import Residue
from Pos import Pos

MAX_TRANSLATION = 10

def cluster_move(main_lattice: Lattice) -> int:
    """ Perform a cluster move on the lattice."""
    cl = Cluster()
    cl.startCounting(main_lattice)
    num_clusters = cl.getTotalClusters()
    
    if num_clusters == 0:
        return -1
        
    # Pick random cluster
    cluster_to_move = random.randint(0, num_clusters - 1)
    chain_list = cl.getChains(cluster_to_move)
    
    # Set rotation
    rotational_dir = random.randint(0, 5)
    
    # Choose random pivot
    x_pivot = random.randint(0, main_lattice.LX - 1)
    y_pivot = random.randint(0, main_lattice.LY - 1)
    z_pivot = random.randint(0, main_lattice.LZ - 1)
    pos_pivot = Pos(x_pivot, y_pivot, z_pivot)
    
    # Store new positions
    new_positions = []
    
    # Check for clashes
    for i, cn in enumerate(chain_list):
        ch = main_lattice.chains[cn]
        chain_positions = []
        for n in range(ch.N):
            res = ch.residues[n]
            new_pos = ch.rotationPosition(res.pos, rotational_dir, pos_pivot)
            chain_positions.append(new_pos)
            if main_lattice.checkForClashesAndTouches(new_pos):
                return -1
        new_positions.append(chain_positions)
    
    # Perform the move
    for i, cn in enumerate(chain_list):
        ch = main_lattice.chains[cn]
        for n in range(ch.N):
            res = ch.residues[n]
            new_spin = ch.rotationSpin(res.spin, rotational_dir)
            main_lattice.emptyLatticePos(res.pos)
            
            res.pos = new_positions[i][n]
            res.spin = new_spin
            main_lattice.setResidue(new_positions[i][n], res)
    
    return 1
