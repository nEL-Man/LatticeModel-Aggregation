import random
from typing import Optional
from GranCan import GranCan
from Lattice import Lattice
from Moves import globalMove, changeStrandCoil, shuffleSpin, localMove
from Cluster import Cluster
from App import App

class MonteCarlo:
    def __init__(self):
        self.main_lattice = Lattice()
        self.grancan = GranCan(self.main_lattice)
        
        self.p_global = 0.1
        self.p_ins_del = 0.01
        self.p_cluster_move = 0.005
        self.p_change_strand_coil = 0.5
        
        self.steps = 10000
        self.total_swaps = 1
    
    def getCext(self) -> int:
        self.main_lattice.stats.printCout()
        print(f"CEXT {App.fn_in} {self.main_lattice.stats.getCext()}")
        return self.main_lattice.stats.getCext()
    
    def MC(self, beta: float, beta_id: int, beta_changed: bool) -> int:
        self.grancan.total_empty_steps = 0
        self.main_lattice.setBetaMoves(beta)
        
        if App.eBetaMu > 0:
            self.grancan.molbox.createNewSet(beta, App.num_configs, beta_changed)
        
        self.main_lattice.energy_map.setBetaID(beta_id)
        
        for _ in range(self.steps):
            if random.random() < self.p_global:
                c = random.choice(self.main_lattice.chains[:self.main_lattice.nChains])
                globalMove(self.main_lattice, c)
            
            if random.random() < self.p_ins_del:
                if random.random() < 0.5:
                    self.grancan.trialInsertChain()
                else:
                    self.grancan.trialDeleteChain()
            
            if random.random() < self.p_change_strand_coil:
                c = random.choice(self.main_lattice.chains[:self.main_lattice.nChains])
                changeStrandCoil(self.main_lattice, c)
            else:
                c = random.choice(self.main_lattice.chains[:self.main_lattice.nChains])
                shuffleSpin(self.main_lattice, c)
            
            c = random.choice(self.main_lattice.chains[:self.main_lattice.nChains])
            localMove(self.main_lattice, c)
        
        if App.clusterMoves and random.random() < self.p_cluster_move:
            self.grancan.clusterMove()
        
        # Check chain constraints
        for c in self.main_lattice.chains[:self.main_lattice.nChains]:
            c.checkAllFwdBkwd()
        
        ci = None
        if self.p_ins_del > 0 or App.clustStats:
            cl = Cluster()
            ci = cl.checkCluster(self.main_lattice, beta)
        
        self.main_lattice.checkStats()
        
        if App.recordMovie:
            self.main_lattice.writePDBMultiChain(App.moviestream, App.moviestep)
            App.moviestep += 1
        
        self.main_lattice.energy_map.writeTableStats(
            self.main_lattice.stats, App.tableStep, 1.0, ci
        )
        
        if self.grancan.total_empty_steps > 0:
            empty_stats = Stats()
            weight = self.grancan.total_empty_steps / float(self.steps)
            ci_empty = ClusterInfo(0, 0, 0)
            self.main_lattice.energy_map.writeTableStats(
                empty_stats, -App.tableStep, weight, ci_empty
            )
        
        App.tableStep += 1
        return self.main_lattice.stats.getEtot()
    
    def getMCStats(self, num_procs: int, my_rank: int):
        cl = Cluster()
        cl.startCounting(self.main_lattice)
        cl.printCout()
        
        # Write output files
        out_pdb = f"outN{my_rank}.pdb"
        self.main_lattice.writePDBMultiChain(out_pdb)
        
        out_stats = f"{App.fn_stats}N{my_rank}.tab"
        self.main_lattice.energy_map.printGnuPlotDataWR(out_stats, num_procs, my_rank)
        
        out_emap = f"{App.fn_eMap}N{my_rank}.tab"
        self.main_lattice.energy_map.printEMap(out_emap, num_procs, my_rank)
        
        if App.hbondStats:
            out_hbond = f"{App.fn_hbond}N{my_rank}.tab"
            self.main_lattice.energy_map.printHBondMap(out_hbond, num_procs, my_rank)
        
        self.main_lattice.stats.printCout()
        print(f"free chains: {self.grancan.getNumFreeChains()}")
        
        self.main_lattice.energy_map.finalizeTableStats() 
