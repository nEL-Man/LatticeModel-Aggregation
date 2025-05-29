from typing import List, Dict, Tuple
from dataclasses import dataclass
import sys
from Lattice import Lattice
from Residue import Residue
from AA import AA
from App import App
import math

MAX_CHAINS = 1000;

@dataclass
class ClusterInfo: 
	clusters: int = 0;
	largestCluster: int = 0; 
	totalChains: int = 0;

class Cluster:
	oldClusterSize = 0;

	def_init_(self):
		self.largestCluster = 0;
		self.totalClusters = 0;
		self.totalVisited = 0;
		self.visited = [-1] * MAX_CHAINS;
		self.clustNum = [-1] * MAX_CHAINS;
		self.clusterCount = [0] * MAX_CHAINS;
		self.1 = None;

	def reset(self):
		self.largestCluster = 0;
		self.totalClusters = 0;
		self.totalVisited = 0
        self.visited = [-1] * MAX_CHAINS
        self.clustNum = [-1] * MAX_CHAINS
        self.clusterCount = [0] * MAX_CHAINS
    
    def startCounting(self, l1: Lattice):
        self.l = l1
        self.reset()
        curCluster = 0
        self.totalChains = l1.nChains
        
        for i in range(l1.nChains - 1, -1, -1):
            if self.visited[i] < 0:
                self.visit(i, curCluster)
                if self.clusterCount[curCluster] > self.clusterCount[self.largestCluster]:
                    self.largestCluster = curCluster
                curCluster += 1
        self.totalClusters = curCluster
    
    def visit(self, chain: int, cluster: int):
        self.clusterCount[cluster] += 1
        self.totalVisited += 1
        self.visited[chain] = cluster
        self.clustNum[chain] = cluster
        
        if self.totalVisited >= self.l.nChains:
            return
            
        for n in range(self.l.chains[chain].N):
            res = self.l.chains[chain].residues[n]
            for k in range(6):
                nb = Pos.local[k] + res.pos
                nb.periodicBoundary()
                neighbour = self.l.getResidue(nb)
                if neighbour is not None and neighbour != self.l.AIR:
                    if not res.masked and not neighbour.masked:
                        nChain = neighbour.chainNum
                        if self.visited[nChain] < 0:
                            self.visit(nChain, cluster)
    
    def checkCluster(self, l1: Lattice, beta: float) -> ClusterInfo:
        self.l = l1
        self.startCounting(l1)
        ci = ClusterInfo()
        
        lowLim = int(0.95 * self.oldClusterSize)
        highLim = int(1.05 * self.oldClusterSize)
        
        if (self.clusterCount[self.largestCluster] < lowLim or 
            self.clusterCount[self.largestCluster] > highLim):
            print(f"Largest cluster {self.largestCluster} size {self.clusterCount[self.largestCluster]}")
            if App.writeCluster:
                self.writePDBClust(self.largestCluster, beta)
            self.oldClusterSize = self.clusterCount[self.largestCluster]
        
        ci.clusters = self.totalClusters
        ci.largestCluster = self.clusterCount[self.largestCluster]
        ci.totalChains = l1.nChains
        return ci
    
    def writePDBClust(self, clust: int, beta: float):
        stT = f"{0.01/beta:.2f}"
        size = self.clusterCount[clust]
        
        # Write full PDB
        outPDB = f"clusT{stT}S{size}all.pdb"
        self.l.writePDBMultiChain(outPDB)
        
        # Write cluster PDB
        outPDBclust = f"clusT{stT}S{size}.pdb"
        with open(outPDBclust, 'w') as outs:
            outs.write("MODEL\n")
            for nc in range(self.l.nChains):
                if self.visited[nc] == clust:
                    ch = self.l.chains[nc]
                    if ch.frozen and ch.locked:
                        outs.write("LIG\n")
                    elif ch.frozen:
                        outs.write("RIG\n")
                        
                    for n in range(ch.N):
                        res = ch.residues[n]
                        # Write atom records...
                        # Similar formatting as in C++ version
                        pass
        
        # Write cluster info file
        outClust = f"clusT{stT}S{size}.out"
        with open(outClust, 'w') as outfile:
            outfile.write(f"beta: {beta}\n")
            outfile.write(f"beta: {stT}\n")
            outfile.write(f"eBetaMu: {App.eBetaMu}\n")
            outfile.write(f"LX LY LZ: {self.l.LX} {self.l.LY} {self.l.LZ}\n")
            outfile.write(self.l.stats.print2string() + "\n")
            
            for cc in range(self.totalClusters):
                outfile.write(f"cluster {cc} contains {self.clusterCount[cc]} chains\n")
    
    def getChains(self, cluster_number: int) -> List[int]:
        return [nc for nc in range(self.totalChains) 
                if cluster_number == self.clustNum[nc]]
    
    def getTotalClusters(self) -> int:
        return self.totalClusters
		
