from typing import Dict, List
import math
from Stats import Stats
from Lattice import Lattice

QMAX = 10000
MAX_PROCS = 100  # Adjust based on your needs

class ClusterInfo:
    def __init__(self):
        self.clusters = 0
        self.largestCluster = 0
        self.totalChains = 0

class EnergyMap:
    def __init__(self, lattice: Lattice):
        self.l = lattice
        self.energy_map = [[0.0 for _ in range(QMAX)] for _ in range(MAX_PROCS)]
        self.hbond_map = [[0.0 for _ in range(QMAX)] for _ in range(MAX_PROCS)]
        self.e_map = [dict() for _ in range(MAX_PROCS)]  # Using dict instead of unordered_map
        self.beta_id = 0
        self.table_stats_stream = None
        self.lowest_efound = 8000
    
    def setBetaID(self, b_id: int):
        self.beta_id = b_id
    
    def mapStats(self, stats: Stats, boltz: float):
        if not App.nativeStats:
            if stats.Cext < 0 or stats.Cext >= QMAX:
                print(f"ERROR too many ({stats.Cext}) external contacts")
                sys.exit(1)
            if math.isnan(boltz):
                print("ERROR boltz nan")
                sys.exit(1)
            self.energy_map[self.beta_id][stats.Cext] += boltz
        else:
            if stats.Nint < 0 or stats.Nint >= QMAX:
                print(f"ERROR too many ({stats.Nint}) native contacts")
                sys.exit(1)
            
            if App.checkBelowNative and self.l.native.isBelowNativeEnergy(stats.Etot):
                self.checkBelowNative(stats, boltz)
            
            if App.cintStats:
                self.energy_map[self.beta_id][stats.Cint] += boltz
            else:
                self.energy_map[self.beta_id][stats.Nint] += boltz
        
        if App.hbondStats:
            self.hbond_map[self.beta_id][stats.Hext] += boltz
        
        self.e_map[self.beta_id][stats.Etot] = self.e_map[self.beta_id].get(stats.Etot, 0) + boltz
    
    def initTableStats(self):
        filename = f"statsTable{App.myRank}.txt"
        self.table_stats_stream = open(filename, 'w')
        headers = [
            "#step", "T", "beta", "weight",
            "Eint", "Eext", "Esol", "Etot",
            "Cint", "Cext", "Ctot"
        ]
        
        if App.nativeStats:
            headers.extend(["Nint", "Next", "Ntot"])
        if App.hbondStats:
            headers.extend(["Hint", "Hext", "Htot"])
        if App.eBetaMu > 0 or App.clustStats:
            headers.extend(["eBM", "clust", "maxCl", "nChains"])
        
        self.table_stats_stream.write("\t".join(headers) + "\n")
    
    def writeTableStats(self, stats: Stats, step: int, weight: float, ci: ClusterInfo):
        beta = self.l.betas[self.beta_id] * 100
        temp = 1 / beta if beta != 0 else float('inf')
        
        data = [
            str(step), f"{temp:.4f}", str(beta), str(weight),
            str(stats.Eint), str(stats.Eext), str(stats.Esol), str(stats.Etot),
            str(stats.Cint), str(stats.Cext), str(stats.Ctot)
        ]
        
        if App.nativeStats:
            data.extend([str(stats.Nint), str(stats.Next), str(stats.Ntot)])
        if App.hbondStats:
            data.extend([str(stats.Hint), str(stats.Hext), str(stats.Htot)])
        if App.eBetaMu > 0 or App.clustStats:
            data.extend([
                str(App.eBetaMu), str(ci.clusters),
                str(ci.largestCluster), str(ci.totalChains)
            ])
        
        self.table_stats_stream.write("\t".join(data) + "\n")
    
    def finalizeTableStats(self):
        if self.table_stats_stream:
            self.table_stats_stream.close()
    
    def checkBelowNative(self, stats: Stats, boltz: float):
        if App.findLowestE:
            if self.lowest_efound > stats.Etot:
                print("energy below native energy")
                print("Mapped stats:")
                stats.printCout()
                print("Lattice stats:")
                self.l.stats.printCout()
                
                if self.l.stats.Etot == stats.Etot:
                    filename = f"lowestEnergyN{App.myRank}.pdb"
                    self.l.writePDBMultiChain(filename)
                    self.lowest_efound = stats.Etot
                    
                    info_filename = f"lowestEnergyN{App.myRank}.txt"
                    with open(info_filename, 'w') as f:
                        f.write(stats.print2string() + "\n")
                        f.write(f"beta {self.l.getBetaMoves()}\n")
        else:
            if not self.l.native.hasNativeStructure(stats.Ntot):
                print("ERROR energy below native energy")
                print("Mapped stats:")
                stats.printCout()
                print("Lattice stats:")
                self.l.stats.printCout()
                
                if self.l.stats.Etot == stats.Etot:
                    self.l.writePDBMultiChain("belowNativeEnergy.pdb")
                    sys.exit(1)
    
    def printGnuPlotDataWR(self, filename: str, num_procs: int, my_rank: int):
        with open(filename, 'w') as f:
            print(f"printing to {filename}")
            for beta_id in range(num_procs):
                f.write(f"BETA\t{100.0 * self.l.betas[beta_id]}\n")
                for q in range(QMAX):
                    if self.energy_map[beta_id][q] > 0:
                        f.write(f"{q}\t{self.energy_map[beta_id][q]}\n")
            print(f"finished printing to {filename}")
    
    def printEMap(self, filename: str, num_procs: int, my_rank: int):
        with open(filename, 'w') as f:
            print(f"printing to {filename}")
            for beta_id in range(num_procs):
                f.write(f"BETA\t{100.0 * self.l.betas[beta_id]}\n")
                for energy, count in self.e_map[beta_id].items():
                    f.write(f"{energy / 100.0}\t{count}\n")
    
    def printHBondMap(self, filename: str, num_procs: int, my_rank: int):
        with open(filename, 'w') as f:
            print(f"printing to {filename}")
            for beta_id in range(num_procs):
                f.write(f"BETA\t{100.0 * self.l.betas[beta_id]}\n")
                for q in range(QMAX):
                    if self.hbond_map[beta_id][q] > 0:
                        f.write(f"{q}\t{self.hbond_map[beta_id][q]}\n")
            print(f"finished printing to {filename}")
