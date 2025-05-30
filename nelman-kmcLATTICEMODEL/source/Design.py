from typing import List, Dict, Optional
import random
import math
import sys
from AA import AA
from Lattice import Lattice
from AADistr import AADistr
from Stats import Stats
from Moves import changeStrandCoil, shuffleSpin

MAXAA = 21

class Design:
    def __init__(self):
        self.main_lattice = Lattice()
        self.main_lattice.setAA(App.fn_aa, App.indx_water)
        self.aa_distr = None
        self.count_aa = [0] * MAXAA
        self.system_var = 0
        self.num_res = 0
        
        if App.AAdistribution:
            self.aa_distr = AADistr()
            self.aa_distr.init(App.fn_AAdistribution, self.main_lattice.aaInt)
    
    def init(self, fn_in: str):
        self.main_lattice.readPDBMultiChain(fn_in)
        if App.AAdistribution:
            self.aa_distr.setSequence(self.main_lattice)
        else:
            self.initCountAA()
        
        # Set random spins
        for cn in range(self.main_lattice.nChains):
            if not self.main_lattice.chains[cn].frozen and App.flipSpins:
                self.main_lattice.chains[cn].setRandomSpins()
    
    def designProcedure(self, beta: float, num_steps: int):
        self.main_lattice.setBetaMoves(beta)
        print(f"energy Beta {beta}")
        
        for i in range(num_steps):
            # Change AA
            cn = random.randint(0, self.main_lattice.nChains - 1)
            if not self.main_lattice.chains[cn].frozen:
                self.changeAA(cn)
                
                if App.flipSpins and random.random() < 0.5:
                    if random.random() < 0.5:
                        # Change state
                        c = self.main_lattice.chains[random.randint(0, self.main_lattice.nChains - 1)]
                        changeStrandCoil(self.main_lattice, c)
                    else:
                        # Change spin
                        c = self.main_lattice.chains[random.randint(0, self.main_lattice.nChains - 1)]
                        shuffleSpin(self.main_lattice, c)
    
    def changeAA(self, chain: int):
        c = self.main_lattice.chains[chain]
        n = random.randint(0, c.N - 1)
        res = c.residues[n]
        old_aa = res.aa
        new_aa = old_aa
        while new_aa == old_aa or new_aa == AA.designWater:
            new_aa = random.randint(0, AA.NUMAA - 1)
        
        # Check distribution constraints
        local_var = 0
        dist = 0
        if App.AAdistribution:
            dist = self.aa_distr.getDistanceSubst(old_aa, new_aa)
            acc = math.exp(-dist / App.designT)
            if random.random() >= acc:
                return
        else:
            local_var = self.count_aa[old_aa] / (self.count_aa[new_aa] + 1.0)
            acc_var = pow(local_var, App.designT)
            if random.random() >= acc_var:
                return
        
        # Calculate energy changes
        old_stats = Stats()
        new_stats = Stats()
        old_stats.localStats(res, res.pos, res.spin, res.state, old_aa, self.main_lattice)
        new_stats.localStats(res, res.pos, res.spin, res.state, new_aa, self.main_lattice)
        
        dE = new_stats.getDeltaE(old_stats)
        boltz = math.exp(-dE * self.main_lattice.getBetaMoves())
        accept = dE <= 0 or random.random() < boltz
        
        if accept:
            res.aa = new_aa
            self.main_lattice.stats = self.main_lattice.stats.delta(new_stats, old_stats)
            
            if App.AAdistribution:
                self.aa_distr.substitute(old_aa, new_aa, dist)
            else:
                self.count_aa[old_aa] -= 1
                self.count_aa[new_aa] += 1
                self.system_var *= local_var
    
    def initCountAA(self):
        self.num_res = 0
        self.count_aa = [0] * AA.NUMAA
        
        for nc in range(self.main_lattice.nChains):
            ch = self.main_lattice.chains[nc]
            for n in range(ch.N):
                res = ch.residues[n]
                self.count_aa[res.aa] += 1
                self.num_res += 1
        
        self.system_var = self.getSystemVar()
    
    def getSystemVar(self) -> float:
        ni = 1.0
        for aa in range(AA.NUMAA):
            if self.count_aa[aa]:
                ni *= math.factorial(self.count_aa[aa])
        return 1.0 / ni
    
    def writeDesignStats(self, filename: str):
        with open(filename, 'w') as f:
            f.write(f"designT:{App.designT}\n")
            f.write(f"energyT:{self.main_lattice.getBetaMoves()}\n")
            
            # Write sequence
            for nc in range(self.main_lattice.nChains):
                ch = self.main_lattice.chains[nc]
                f.write(f"chain no:{nc}\n")
                for n in range(ch.N):
                    res = ch.residues[n]
                    f.write(f"{AA.int2aa[res.aa]} ")
                f.write("\n")
            
            # Write stats
            f.write(self.main_lattice.stats.print2string())
            
            # Write distribution info
            if App.AAdistribution:
                f.write(self.aa_distr.toString())
            else:
                f.write(f"systemsVar:{self.system_var}\n")
                f.write("AA count:\n")
                for aa in range(AA.NUMAA):
                    f.write(f"{AA.int2aa[aa]} {self.count_aa[aa]}, ")
                f.write("\n")
