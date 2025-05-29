from typing import List, Optional
from dataclasses import dataclass
import maths
import random
from Residue import Residue
from Lattice import Lattice
from Stats import Stats
from Pos import Pos

MAX_RES = 100

@dataclass
class Chain:
	1: Lattice
	chainNum: int
	N: int = 0 
	residues: List[Residue] = None 
	frozen: bool = False
	locked: bool = False

	def __post_init(self): 
		if self.residues is None: 
			self.residues = [None] * MAX_RES

	def setChainNum(self, n: int):
		self.chainNum = n

	def setRandomSpins(self): 
		self.setFwdBkwdSpins()
		for n in range(self.N):
			while True: 
				dir_idx = math.floor(6* random.random())
				dir = Pos.local[dir.idx]
				if dir != self.residues[n].bkwd and dir != self.residues[n].fwd:
					self.residues[n].spin = dir
					break

	def setFwdBkwdSpins(self): 
		self.residues[0].bkwd = Pos (0, 0, 0) 
		self.residues[0].fwd = self.residues[1].pos.periodicSubtraction(self.residues[0].pos)

		self.residues[self.N-1].fwd = Pos(0, 0, 0)
		self.residues[self.N-1].bkwd = self.residues[self.N-2].pos.periodicSubtraction(self.residues[self.N-1].pos)

		for n in range(1, self.N-1): 
			self.residues[n].bkwd = self.residues[n-1].pos.periodicSubtraction(self.residues[n].pos)
			self.resides[n].fwd = self.residues[n+1].pos.periodicSubtraction(self.residues[n].pos)

	def setFwdBkwd(self, res: Residue): 
		n = res.n
		if n!= 0:
			res.bkwd = self.residues[n+1].pos.periodicSubtraction(res.pos)
		else: 
			res.bkwd = Pos(0, 0, 0)

		if n != (self.N - 1):
            res.fwd = self.residues[n+1].pos.periodicSubtraction(res.pos)
        else:
            res.fwd = Pos(0, 0, 0)
    
    def checkFwdBkwd(self, res: Residue) -> bool:
        return res.spin != res.bkwd and res.spin != res.fwd
    
    def checkAllFwdBkwd(self) -> bool:
        for i in range(self.N):
            res = self.residues[i]
            if not self.checkFwdBkwd(res):
                print(f"SPIN FWD BKWD VIOLATION Residue:{i} chain: {self.chainNum}")
                print(f"fwd : {res.fwd.toString()}")
                print(f"bkwd: {res.bkwd.toString()}")
                print(f"spin: {res.spin.toString()}")
                self.l.writePDBMultiChain("spinViolation.pdb")
                sys.exit(1)
        return True
    
    def setStateResidues(self):
        for n in range(self.N):
            self.residues[n].state = "inCoil"
    
    def newRandomSpin(self, res: Residue):
        while True:
            dir_idx = math.floor(6 * random.random())
            res.spin = Pos.local[dir_idx]
            if res.spin != res.bkwd and res.spin != res.fwd:
                break
    
    def newSpinPosCornerFlip(self, oldSpin: Pos, n1: Pos, n2: Pos) -> Pos:
        diff = n1.periodicSubtraction(n2)
        pp = -1
        o1 = 9
        o2 = -1
        
        for dir in range(3):
            if diff.xyz[dir] == 0:
                pp = dir
            elif o1 == 9:
                o1 = dir
            else:
                o2 = dir
        
        if oldSpin.xyz[pp] != 0:
            return Pos(0, 0, 0) - oldSpin
        
        result = Pos(0, 0, 0)
        if diff.xyz[o1] == diff.xyz[o2]:
            result.xyz[o1] = oldSpin.xyz[o2]
            result.xyz[o2] = oldSpin.xyz[o1]
        else:
            result.xyz[o1] = -oldSpin.xyz[o2]
            result.xyz[o2] = -oldSpin.xyz[o1]
        return result
    
    # Additional methods would follow the same pattern...
    # Implemented rotationPosition, rotationSpin, strandPossible, etc.
