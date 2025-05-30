from typing import List, Optional
import random
from Lattice import Lattice
from Chain import Chain
from Config import Config
from Moves import localMove, rotatePoint, changeStrandCoil, shuffleSpin

MAX_CONFIGS = 1000
EQUIL_STEPS = 1000

class MolBox:
    def __init__(self, main_l: Lattice):
        self.current_config = 0
        self.num_configs = 1000
        self.box = Lattice()
        self.main_lattice = main_l
        self.config_set = [None] * MAX_CONFIGS
        self.sequence = [0] * MAX_RES
        self.molecule = None
    
    def setAA(self, fn: str, indx: int):
        self.box.setAA(fn, indx)
    
    def getBeta(self) -> float:
        return self.box.getBetaMoves()
    
    def setMolecule(self, chain: Chain):
        self.molecule = Chain(self.box, 0)
        self.box.chains[0] = self.molecule
        self.molecule.frozen = False
        self.molecule.locked = False
        self.molecule.N = chain.N
        
        for n in range(chain.N):
            res = Residue()
            res.pos = chain.residues[n].pos.copy()
            res.spin = chain.residues[n].spin.copy()
            res.state = chain.residues[n].state
            res.aa = chain.residues[n].aa
            self.sequence[n] = res.aa
            res.n = chain.residues[n].n
            res.chainNum = 0
            self.molecule.residues[n] = res
            self.box.r[res.pos.x][res.pos.y][res.pos.z] = res
        
        self.molecule.setFwdBkwdSpins()
        self.box.nChains = 1
        self.box.energyMap = None
    
    def createNewSet(self, b: float, num_configs: int, beta_changed: bool):
        self.box.setBetaMoves(b)
        
        if beta_changed:
            self.equilibrate(EQUIL_STEPS)
        
        self.current_config = 0
        self.num_configs = min(num_configs, MAX_CONFIGS)
        
        for nc in range(self.num_configs):
            self.equilibrate(random.randint(0, 99))
            self.recordConfig(nc)
    
    def getNewConfig(self) -> Optional[Config]:
        self.current_config += 1
        if self.current_config >= self.num_configs:
            self.createNewSet(self.box.getBetaMoves(), self.num_configs, False)
        return self.config_set[random.randint(0, self.num_configs - 1)]
    
    def equilibrate(self, n_steps: int):
        for _ in range(n_steps):
            localMove(self.box, self.molecule)
            if random.random() < 0.05:
                rotatePoint(self.box, self.molecule)
            if random.random() < 0.5:
                changeStrandCoil(self.box, self.molecule)
            else:
                shuffleSpin(self.box, self.molecule)
    
    def recordConfig(self, insert_at: int):
        if self.config_set[insert_at] is not None:
            del self.config_set[insert_at]
        
        cc = Config()
        for n in range(self.molecule.N):
            cc.positions[n] = self.molecule.residues[n].pos.copy()
            cc.spins[n] = self.molecule.residues[n].spin.copy()
            cc.states[n] = self.molecule.residues[n].state
            cc.nTotal = self.molecule.N
        
        self.config_set[insert_at] = cc
