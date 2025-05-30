from enum import Enum
from typing import List, Optional, Any
import sys

class State(Enum):
    inCoil = 0
    inStrand = 1

class Pos:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z
    
    @staticmethod
    def orthogonal(p1, p2):
        return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z == 0

class Residue:
    def __init__(self):
        self.aa = 0
        self.pos = Pos()
        self.spin = Pos()
        self.bkwd = Pos()
        self.fwd = Pos()
        self.n = 0
        self.chainNum = 0
        self.state = State.inCoil
        self.masked = False

class Chain:
    def __init__(self):
        self.N = 0
        self.residues = []
        self.chainID = ""

class Lattice:
    def __init__(self):
        self.nChains = 0
        self.chains = []
        self.AIR = None
        self.stats = None
        self.native = None
        self.aaInt = None
    
    def getResidue(self, pos):
        for chain in self.chains:
            for res in chain.residues:
                if res.pos == pos:
                    return res
        return None

class App:
    nativeStats = False
    hbondStats = False
    BBSolvE = 0
    HBondE = 0
    StericE = 0

# Global variables
local = [Pos() for _ in range(6)]
cross = [Pos() for _ in range(12)]

class Stats:
    if hasattr(sys, 'CROSSNBS'):
        loopUpCrossNbs = [[[[Pos() for _ in range(4)] for _ in range(3)] for _ in range(3)] for _ in range(3)]
    
    def __init__(self):
        self.Etot = 0
        self.Ctot = 0
        self.Cext = 0
        self.Eint = 0
        self.Cint = 0
        self.Eext = 0
        self.Esol = 0
        self.Nint = 0
        self.Next = 0
        self.Ntot = 0
        self.Hint = 0
        self.Hext = 0
        self.Htot = 0
    
    def clean(self):
        self.Esol = 0
        self.Eext = 0
        self.Eint = 0
        self.Etot = 0
        self.Cint = 0
        self.Cext = 0
        self.Ctot = 0
        if App.nativeStats:
            self.Nint = 0
            self.Next = 0
            self.Ntot = 0
        if App.hbondStats:
            self.Hint = 0
            self.Hext = 0
            self.Htot = 0
    
    def __add__(self, other):
        result = Stats()
        result.Cint = self.Cint + other.Cint
        result.Cext = self.Cext + other.Cext
        result.Ctot = self.Ctot + other.Ctot
        result.Eint = self.Eint + other.Eint
        result.Eext = self.Eext + other.Eext
        result.Etot = self.Etot + other.Etot
        result.Esol = self.Esol + other.Esol
        
        if App.nativeStats:
            result.Nint = self.Nint + other.Nint
            result.Next = self.Next + other.Next
            result.Ntot = self.Ntot + other.Ntot
        
        if App.hbondStats:
            result.Hint = self.Hint + other.Hint
            result.Hext = self.Hext + other.Hext
            result.Htot = self.Htot + other.Htot
        
        return result
    
    def __iadd__(self, other):
        self.Cint += other.Cint
        self.Cext += other.Cext
        self.Ctot += other.Ctot
        self.Eint += other.Eint
        self.Eext += other.Eext
        self.Etot += other.Etot
        self.Esol += other.Esol
        
        if App.nativeStats:
            self.Nint += other.Nint
            self.Next += other.Next
            self.Ntot += other.Ntot
        
        if App.hbondStats:
            self.Hint += other.Hint
            self.Hext += other.Hext
            self.Htot += other.Htot
        
        return self
    
    def __isub__(self, other):
        self.Cint -= other.Cint
        self.Cext -= other.Cext
        self.Ctot -= other.Ctot
        self.Eint -= other.Eint
        self.Eext -= other.Eext
        self.Etot -= other.Etot
        self.Esol -= other.Esol
        
        if App.nativeStats:
            self.Nint -= other.Nint
            self.Next -= other.Next
            self.Ntot -= other.Ntot
        
        if App.hbondStats:
            self.Hint -= other.Hint
            self.Hext -= other.Hext
            self.Htot -= other.Htot
        
        return self
    
    def __eq__(self, other):
        return not (self != other)
    
    def __ne__(self, other):
        ans = (self.Cint != other.Cint or
               self.Cext != other.Cext or
               self.Ctot != other.Ctot or
               self.Eint != other.Eint or
               self.Eext != other.Eext or
               self.Esol != other.Esol or
               self.Etot != other.Etot)
        
        if App.nativeStats:
            ans = (ans or
                   self.Nint != other.Nint or
                   self.Next != other.Next or
                   self.Ntot != other.Ntot)
        
        if App.hbondStats:
            ans = (ans or
                   self.Hint != other.Hint or
                   self.Hext != other.Hext or
                   self.Htot != other.Htot)
        
        return ans
    
    def delta(self, Snew, Sold):
        out = Stats()
        out.Eext = self.Eext + Snew.Eext - Sold.Eext
        out.Eint = self.Eint + Snew.Eint - Sold.Eint
        out.Etot = self.Etot + Snew.Etot - Sold.Etot
        out.Esol = self.Esol + Snew.Esol - Sold.Esol

        out.Cext = self.Cext + Snew.Cext - Sold.Cext
        out.Cint = self.Cint + Snew.Cint - Sold.Cint
        out.Ctot = self.Ctot + Snew.Ctot - Sold.Ctot

        if App.nativeStats:	 
            out.Next = self.Next + Snew.Next - Sold.Next
            out.Nint = self.Nint + Snew.Nint - Sold.Nint
            out.Ntot = self.Ntot + Snew.Ntot - Sold.Ntot
        
        if App.hbondStats:	 
            out.Hext = self.Hext + Snew.Hext - Sold.Hext
            out.Hint = self.Hint + Snew.Hint - Sold.Hint
            out.Htot = self.Htot + Snew.Htot - Sold.Htot

        return out
    
    def getDeltaE(self, old):
        return self.Etot - old.Etot
    
    def getCext(self):
        return self.Cext
    
    def getEtot(self):
        return self.Etot
    
    @classmethod
    def setLoopUpCrossNbs(cls):
        cnt = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
        
        for k in range(6):
            spin = local[k]
            for i in range(12):
                if ((spin.x != 0 and cross[i].x == spin.x) or 
                    (spin.y != 0 and cross[i].y == spin.y) or
                    (spin.z != 0 and cross[i].z == spin.z)):
                    indx = cnt[spin.x+1][spin.y+1][spin.z+1]
                    cls.loopUpCrossNbs[spin.x+1][spin.y+1][spin.z+1][indx] = cross[i]
                    cnt[spin.x+1][spin.y+1][spin.z+1] += 1
    
    @classmethod
    def getCrossNbsSpinDir(cls, spin):
        return cls.loopUpCrossNbs[spin.x+1][spin.y+1][spin.z+1]
    
    def getLatticeStats(self, l):
        self.Eext = 0
        self.Eint = 0
        self.Cint = 0
        self.Cext = 0
        self.Ctot = 0
        self.Etot = 0
        self.Esol = 0
        self.Nint = 0
        self.Next = 0
        self.Ntot = 0
        self.Hint = 0
        self.Hext = 0
        self.Htot = 0

        for cn in range(l.nChains):
            c = l.chains[cn]
            for n in range(c.N):
                res = c.residues[n]
                tmp = Stats()
                tmp.localStats(res, res.pos, res.spin, res.state, l)
                self += tmp
        
        self.Eext //= 2
        self.Eint //= 2
        self.Cint //= 2
        self.Cext //= 2

        self.Etot = self.Eext + self.Eint + self.Esol
        self.Ctot = self.Cint + self.Cext

        self.Nint //= 2
        self.Next //= 2
        self.Ntot = self.Nint + self.Next
        self.Hint //= 2
        self.Hext //= 2
        self.Htot = self.Hint + self.Hext
    
    def get_Eint_Cint(self, c, l):
        self.Eint = 0
        self.Cint = 0
        self.Esol = 0
        for n in range(c.N):
            res = c.residues[n]
            tmp = Stats()
            tmp.localStats(res, res.pos, res.spin, res.state, l)
            self += tmp
        
        self.Eint //= 2
        self.Cint //= 2
        self.Etot = self.Eint + self.Esol
        self.Ctot = self.Cint
        
        # Set all others to zero
        self.Eext = 0
        self.Cext = 0

        self.Nint //= 2  
        self.Next = 0
        self.Ntot = self.Nint

        self.Hint //= 2  
        self.Hext = 0
        self.Htot = self.Hint
    
    def printCout(self):
        print(f"Eint {self.Eint}", end=' ')
        print(f"Eext {self.Eext}", end=' ')
        print(f"Esol {self.Esol}", end=' ')
        print(f"Etot {self.Etot}")

        print(f"Cint {self.Cint}", end=' ')
        print(f"Cext {self.Cext}", end=' ')
        print(f"Ctot {self.Ctot}")

        if App.nativeStats:	 
            print(f"Nint {self.Nint}", end=' ')
            print(f"Next {self.Next}", end=' ')
            print(f"Ntot {self.Ntot}")
        
        if App.hbondStats: 
            print(f"Hint {self.Hint}", end=' ')
            print(f"Hext {self.Hext}", end=' ')
            print(f"Htot {self.Htot}")
    
    def print2string(self):
        output = []
        output.append(f"Eint {self.Eint}")
        output.append(f"Eext {self.Eext}")
        output.append(f"Esol {self.Esol}")
        output.append(f"Etot {self.Etot}")

        output.append(f"Cint {self.Cint}")
        output.append(f"Cext {self.Cext}")
        output.append(f"Ctot {self.Ctot}")
        
        if App.nativeStats:	  
            output.append(f"Nint {self.Nint}")
            output.append(f"Next {self.Next}")
            output.append(f"Ntot {self.Ntot}")

        if App.hbondStats:   
            output.append(f"Hint {self.Hint}")
            output.append(f"Hext {self.Hext}")
            output.append(f"Htot {self.Htot}")

        return '\n'.join(output)
