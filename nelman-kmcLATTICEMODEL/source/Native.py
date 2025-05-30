from typing import List, Dict, Tuple
import numpy as np

class Native:
    def __init__(self, lattice):
        print("creating new native")
        self.contacts = np.zeros(
            (lattice.nChains, lattice.nChains, 
             max(c.N for c in lattice.chains[:lattice.nChains]),
             max(c.N for c in lattice.chains[:lattice.nChains])),
            dtype=bool
        )
    
    def setNative(self, lattice):
        print("setting native")
        total_c = 0
        
        for nc in range(lattice.nChains):
            for n in range(lattice.chains[nc].N):
                res = lattice.chains[nc].residues[n]
                
                for k in range(6):
                    pos_nb = Pos.local[k] + res.pos
                    pos_nb.periodicBoundary()
                    res_nb = lattice.getResidue(pos_nb)
                    
                    if res_nb is not None:
                        if res_nb.chainNum != res.chainNum or abs(res_nb.n - res.n) != 1:
                            self.contacts[res.chainNum][res_nb.chainNum][res.n][res_nb.n] = True
                            self.contacts[res_nb.chainNum][res.chainNum][res_nb.n][res.n] = True
                            total_c += 1
        
        self.tot_c_nat = total_c // 2
        print(f"TOTAL NATIVE CONTACTS: {self.tot_c_nat}")
        self.native_energy = lattice.stats.getEtot()
        print(f"NATIVE INTERNAL ENERGY: {self.native_energy}")
    
    def isNativeContact(self, chain_a: int, chain_b: int, res_a: int, res_b: int) -> bool:
        return self.contacts[chain_a][chain_b][res_a][res_b]
    
    def isBelowNativeEnergy(self, e_tot: int) -> bool:
        return e_tot < self.native_energy
    
    def getNativeEnergy(self) -> int:
        return self.native_energy
    
    def hasNativeStructure(self, n_tot: int) -> bool:
        return n_tot == self.tot_c_nat
