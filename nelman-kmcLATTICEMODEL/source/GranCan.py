from typing import List, Optional
import random
import math
from Lattice import Lattice
from MolBox import MolBox
from Cluster import Cluster
from Pos import Pos
from Config import Config
from Residue import Residue

MAX_CHAINS = 50000
MAX_TRANSLATION = 10

class GranCan:
    def __init__(self, main_lattice: Lattice):
        self.main_lattice = main_lattice
        self.molbox = MolBox(main_lattice)
        self.volume = main_lattice.LX * main_lattice.LY * main_lattice.LZ
        self.p_global = 0.0
        self.p_ins_del = 0.0
        self.p_change_strand_coil = 0.0
        self.p_cluster_move = 0.0
        self.steps = 0
        self.total_swaps = 0
        self.total_empty_steps = 0.0
    
    def trialInsertChain(self) -> int:
        # Find free chains
        free_chains = [c for c in range(self.main_lattice.nChains) 
                      if not self.main_lattice.chains[c].hasCext()]
        num_free = len(free_chains)
        
        # Calculate acceptance ratio
        prefactor = 1.0 / math.pow(self.main_lattice.getRealBeta(), 1.5) if App.LAMBDA_COR else 1.0
        acc_ratio = prefactor * App.eBetaMu * self.volume / (num_free + 1.0)
        
        if random.random() < acc_ratio:
            config = self.molbox.getNewConfig()
            if config is None:
                return -1
                
            # Add random translation
            add_p = Pos(
                random.randint(0, self.main_lattice.LX - 1),
                random.randint(0, self.main_lattice.LY - 1),
                random.randint(0, self.main_lattice.LZ - 1)
            )
            
            # Check for clashes
            for n in range(config.nTotal):
                pos_n = config.positions[n] + add_p
                pos_n.periodicBoundary()
                
                if self.main_lattice.getResidue(pos_n) is not None:
                    return -1
                    
                for k in range(6):
                    pos_nb = Pos.local[k] + pos_n
                    pos_nb.periodicBoundary()
                    if self.main_lattice.getResidue(pos_nb) is not None:
                        return -1
            
            # Insert the chain
            self.main_lattice.insertChain(config, self.molbox.sequence, self.main_lattice.getBetaMoves())
            return 1
            
        return 0
    
    def trialDeleteChain(self) -> int:
        free_chains = [c for c in range(self.main_lattice.nChains) 
                      if not self.main_lattice.chains[c].hasCext()]
        num_free = len(free_chains)
        
        if num_free < 1:
            return -2
        if self.main_lattice.nChains <= 1:
            return -1
            
        # Calculate acceptance ratio
        prefactor = math.pow(self.main_lattice.getRealBeta(), 1.5) if App.LAMBDA_COR else 1.0
        acc_ratio = prefactor * num_free / (App.eBetaMu * self.volume)
        
        if random.random() < acc_ratio:
            chain = random.choice(free_chains)
            self.main_lattice.deleteChain(chain)
            return 1
            
        return 0
    
    def getNumFreeChains(self) -> int:
        return sum(1 for c in range(self.main_lattice.nChains) 
                 if not self.main_lattice.chains[c].hasCext())
    
    def clusterMove(self) -> int:
        # Similar to cluster_move function above
        cl = Cluster()
        cl.startCounting(self.main_lattice)
        num_clusters = cl.getTotalClusters()
        
        if num_clusters == 0:
            return -1
            
        cluster_to_move = random.randint(0, num_clusters - 1)
        chain_list = cl.getChains(cluster_to_move)
        
        rotational_dir = random.randint(0, 5)
        pivot = Pos(
            random.randint(0, self.main_lattice.LX - 1),
            random.randint(0, self.main_lattice.LY - 1),
            random.randint(0, self.main_lattice.LZ - 1)
        )
        
        # Check for clashes
        new_positions = []
        for cn in chain_list:
            ch = self.main_lattice.chains[cn]
            chain_pos = []
            for n in range(ch.N):
                res = ch.residues[n]
                new_pos = ch.rotationPosition(res.pos, rotational_dir, pivot)
                chain_pos.append(new_pos)
                if self.main_lattice.checkForClashesAndTouches(new_pos):
                    return -1
            new_positions.append(chain_pos)
        
        # Perform the move
        for i, cn in enumerate(chain_list):
            ch = self.main_lattice.chains[cn]
            for n in range(ch.N):
                res = ch.residues[n]
                new_spin = ch.rotationSpin(res.spin, rotational_dir)
                self.main_lattice.emptyLatticePos(res.pos)
                res.pos = new_positions[i][n]
                res.spin = new_spin
                self.main_lattice.setResidue(new_positions[i][n], res)
            ch.setFwdBkwdSpins()
        
        return 1
