import random
import math
from typing import List, Tuple, Optional
from enum import Enum

# Constants
MAX_RES = 1000  # Maximum number of residues
XDIR = 0
YDIR = 1
ZDIR = 2

# Assuming these classes are defined elsewhere
class Pos:
    pass

class Residue:
    pass

class Chain:
    pass

class Lattice:
    pass

class Stats:
    pass

class State(Enum):
    inCoil = 0
    inStrand = 1

# Local directions (assuming these are defined elsewhere)
local = [Pos() for _ in range(6)]

class Moves:
    # Class variables to store changes in stats for speedup
    oldLocalStats = Stats()
    newLocalStats = Stats()
    newLatticeStats = Stats()

    @staticmethod
    def shuffleSpin(l: Lattice, c: Chain) -> int:
        n = random.randint(0, c.N - 1)
        res = c.residues[n]
        
        if res.state == State.inStrand:
            l.sampleCurrentStats()
            return -1
            
        newSpin = local[random.randint(0, 5)]
        if newSpin == res.bkwd or newSpin == res.fwd:
            l.sampleCurrentStats()
            return -1
            
        # Energy acceptance + sampling
        Moves.oldLocalStats.clean()
        Moves.newLocalStats.clean()
        
        Moves.oldLocalStats.localStats(res, res.pos, res.spin, res.state, l)
        Moves.newLocalStats.localStats(res, res.pos, newSpin, res.state, l)
        
        Moves.newLatticeStats.clean()
        Moves.newLatticeStats = l.stats.delta(Moves.newLocalStats, Moves.oldLocalStats)
        
        dE = Moves.newLocalStats.getDeltaE(Moves.oldLocalStats)
        boltz = math.exp(-dE * l.getBetaMoves())
        bb = boltz / (boltz + 1.0)
        
        if l.energyMap is not None:
            l.energyMap.mapStats(l.stats, 1.0 - bb)
            l.energyMap.mapStats(Moves.newLatticeStats, bb)
        
        if dE <= 0 or random.random() < boltz:
            # Accept move
            res.spin = newSpin
            l.stats = Moves.newLatticeStats
            return 1
        else:
            return 0

    @staticmethod
    def changeStrandCoil(l: Lattice, c: Chain) -> int:
        n = random.randint(0, c.N - 1)
        res = c.residues[n]
        
        if res.state == State.inCoil:
            if c.strandPossible(res):
                newState = State.inStrand
            else:
                l.sampleCurrentStats()
                return -1
        else:
            newState = State.inCoil
            
        Moves.oldLocalStats.clean()
        Moves.newLocalStats.clean()
        
        Moves.oldLocalStats.localStats(res, res.pos, res.spin, res.state, l)
        Moves.newLocalStats.localStats(res, res.pos, res.spin, newState, l)
        
        Moves.newLatticeStats.clean()
        Moves.newLatticeStats = l.stats.delta(Moves.newLocalStats, Moves.oldLocalStats)
        
        dE = Moves.newLocalStats.getDeltaE(Moves.oldLocalStats)
        boltz = math.exp(-dE * l.getBetaMoves())
        bb = boltz / (boltz + 1.0)
        
        if l.energyMap is not None:
            l.energyMap.mapStats(l.stats, 1.0 - bb)
            l.energyMap.mapStats(Moves.newLatticeStats, bb)
        
        if dE <= 0 or random.random() < boltz:
            # Accept move
            res.state = newState
            l.stats = Moves.newLatticeStats
            return 1
        else:
            return 0

    @staticmethod
    def localMove(l: Lattice, c: Chain) -> int:
        posNew = Pos()
        newSpin = Pos()
        px, py, pz = -1, -1, -1
        
        n = random.randint(0, c.N - 1)
        
        # Check if move is allowed
        if (c.residues[n].state == State.inStrand or
            (n > 0 and c.residues[n-1].state == State.inStrand) or
            (n < c.N-1 and c.residues[n+1].state == State.inStrand)):
            l.sampleCurrentStats()
            return -1
            
        res = c.residues[n]
        
        if n == 0 or n == c.N - 1:
            # End of chain - choose one of 6 directions
            nn = 1 if n == 0 else c.N - 2
            rr = random.randint(0, 5)
            posNew = local[rr] + c.residues[nn].pos
            posNew.periodicBoundary()
            
            if res.pos != posNew:
                newSpin = c.newSpinEndMove(res.pos, posNew, c.residues[nn].pos, res.spin)
            else:
                newSpin = res.spin
        else:
            # Try corner flip
            px = (res.pos.x == c.residues[n-1].pos.x == c.residues[n+1].pos.x)
            py = (res.pos.y == c.residues[n-1].pos.y == c.residues[n+1].pos.y)
            pz = (res.pos.z == c.residues[n-1].pos.z == c.residues[n+1].pos.z)
            
            if px + py + pz == 2:
                # On a straight line
                l.sampleCurrentStats()
                return -1
            elif px + py + pz == 1:
                if px:
                    posNew.x = res.pos.x
                    posNew.y = c.residues[n-1].pos.y + c.residues[n+1].pos.y - res.pos.y
                    posNew.z = c.residues[n-1].pos.z + c.residues[n+1].pos.z - res.pos.z
                elif py:
                    posNew.y = res.pos.y
                    posNew.x = c.residues[n-1].pos.x + c.residues[n+1].pos.x - res.pos.x
                    posNew.z = c.residues[n-1].pos.z + c.residues[n+1].pos.z - res.pos.z
                else:
                    posNew.z = res.pos.z
                    posNew.x = c.residues[n-1].pos.x + c.residues[n+1].pos.x - res.pos.x
                    posNew.y = c.residues[n-1].pos.y + c.residues[n+1].pos.y - res.pos.y
                
                if res.pos != posNew:
                    newSpin = c.newSpinPosCornerFlip(res.spin, c.residues[n-1].pos, c.residues[n+1].pos)
                else:
                    newSpin = res.spin
            else:
                print(f"Break in chain: {c.chainNum}")
                print(f"Trying to move: {n} with coords:")
                print(f"{res.pos.x} {res.pos.y} {res.pos.z}")
                print(f"{c.residues[n-1].pos.x} {c.residues[n-1].pos.y} {c.residues[n-1].pos.z}")
                print(f"{c.residues[n+1].pos.x} {c.residues[n+1].pos.y} {c.residues[n+1].pos.z}")
                l.writePDBMultiChain("chainBreak.pdb")
                sys.exit(1)
        
        # Check for steric clashes
        res_at_new_pos = l.getResidue(posNew)
        if res_at_new_pos is not None:
            if (res_at_new_pos.n == n-2 and res_at_new_pos.chainNum == c.chainNum):
                return Moves.crankshaft(l, c, n-1, n, px, py, pz)
            elif (res_at_new_pos.n == n+2 and res_at_new_pos.chainNum == c.chainNum):
                return Moves.crankshaft(l, c, n, n+1, px, py, pz)
            else:
                l.sampleCurrentStats()
                return -1
        
        # Check spin directions
        newFwd, newBkwd = Pos(), Pos()
        newFwdPrev, newBkwdNxt = Pos(), Pos()
        
        if n > 0:
            newBkwd.periodicSubtraction(c.residues[n-1].pos, posNew)
            newFwdPrev = Pos(0,0,0) - newBkwd
            if c.residues[n-1].spin == newFwdPrev:
                l.sampleCurrentStats()
                return -1
                
        if n < c.N - 1:
            newFwd.periodicSubtraction(c.residues[n+1].pos, posNew)
            newBkwdNxt = Pos(0,0,0) - newFwd
            if c.residues[n+1].spin == newBkwdNxt:
                l.sampleCurrentStats()
                return -1
        
        # Update lattice
        l.emptyLatticePos(res.pos)
        
        # Calculate old and new stats
        Moves.oldLocalStats.clean()
        Moves.newLocalStats.clean()
        
        Moves.oldLocalStats.localStats(res, res.pos, res.spin, res.state, l)
        Moves.oldLocalStats.solventStats(posNew, l)
        
        Moves.newLocalStats.localStats(res, posNew, newSpin, newFwd, newBkwd, res.state, res.aa, l)
        Moves.newLocalStats.solventStats(res.pos, l)
        
        Moves.newLatticeStats.clean()
        Moves.newLatticeStats = l.stats.delta(Moves.newLocalStats, Moves.oldLocalStats)
        
        dE = Moves.newLocalStats.getDeltaE(Moves.oldLocalStats)
        boltz = math.exp(-dE * l.getBetaMoves())
        bb = boltz / (boltz + 1.0)
        
        if l.energyMap is not None:
            l.energyMap.mapStats(l.stats, 1.0 - bb)
            l.energyMap.mapStats(Moves.newLatticeStats, bb)
        
        if dE <= 0 or random.random() < boltz:
            # Accept move
            if n > 0:
                res.bkwd = newBkwd
                c.residues[n-1].fwd = newFwdPrev
            if n < c.N - 1:
                res.fwd = newFwd
                c.residues[n+1].bkwd = newBkwdNxt
            
            res.spin = newSpin
            res.pos = posNew
            l.setResidue(posNew, res)
            l.stats = Moves.newLatticeStats
            
            return 1 if dE == 0 else 2
        else:
            l.setResidue(res.pos, res)
            return 0

    @staticmethod
    def crankshaft(l: Lattice, c: Chain, n1: int, n2: int, px: int, py: int, pz: int) -> int:
        resn1 = c.residues[n1]
        resn2 = c.residues[n2]
        
        if resn1.state == State.inStrand or resn2.state == State.inStrand:
            l.sampleCurrentStats()
            return -1
            
        if n1 < n2:
            n0 = n1 - 1
            n3 = n1 + 2
        else:
            n0 = n1 + 1
            n3 = n1 - 2
            n0, n3 = n3, n0
            n1, n2 = n2, n1
            
        resn0 = c.residues[n0]
        resn3 = c.residues[n3]
        
        if (n0 > 0 and resn0.state == State.inStrand) or (n3 < c.N-1 and resn3.state == State.inStrand):
            l.sampleCurrentStats()
            return -1
            
        if n0 < 0 or n3 < 0 or n0 >= c.N or n3 >= c.N:
            l.sampleCurrentStats()
            return -3
            
        # Determine rotation directions
        dirOrthogonal, dirStalk, dirLegs = 0, 0, 0
        if px:
            dirOrthogonal = XDIR
            if (resn1.pos.y - resn2.pos.y) == 0:
                dirStalk = YDIR
                dirLegs = ZDIR
            else:
                dirStalk = ZDIR
                dirLegs = YDIR
        elif py:
            dirOrthogonal = YDIR
            if (resn1.pos.x - resn2.pos.x) == 0:
                dirStalk = XDIR
                dirLegs = ZDIR
            else:
                dirStalk = ZDIR
                dirLegs = XDIR
        else:
            dirOrthogonal = ZDIR
            if (resn1.pos.x - resn2.pos.x) == 0:
                dirStalk = XDIR
                dirLegs = YDIR
            else:
                dirStalk = XDIR
                dirLegs = YDIR
                
        dirFlip = 1 if random.random() < 0.5 else -1
        
        # Calculate new positions
        old0 = Pos(resn0.pos.x, resn0.pos.y, resn0.pos.z)
        old3 = Pos(resn3.pos.x, resn3.pos.y, resn3.pos.z)
        
        newPos1 = Pos()
        newPos2 = Pos()
        
        newPos1[dirOrthogonal] = old0[dirOrthogonal] + dirFlip
        newPos2[dirOrthogonal] = old3[dirOrthogonal] + dirFlip
        newPos1[dirStalk] = old0[dirStalk]
        newPos2[dirStalk] = old3[dirStalk]
        newPos1[dirLegs] = old0[dirLegs]
        newPos2[dirLegs] = old3[dirLegs]
        
        newPos1.periodicBoundary()
        newPos2.periodicBoundary()
        
        # Check for collisions
        if l.getResidue(newPos1) is not None or l.getResidue(newPos2) is not None:
            l.sampleCurrentStats()
            return -1
            
        # Check spin directions
        newBkwdn1, newFwdn2 = Pos(), Pos()
        newBkwdn1.periodicSubtraction(resn0.pos, newPos1)
        newFwdn2.periodicSubtraction(resn3.pos, newPos2)
        newFwdn0 = Pos(0,0,0) - newBkwdn1
        newBkwdn3 = Pos(0,0,0) - newFwdn2
        
        if resn0.spin == newFwdn0 or resn3.spin == newBkwdn3:
            l.sampleCurrentStats()
            return -1
            
        newSpin1 = c.newSpinEndMove(resn1.pos, newPos1, resn0.pos, resn1.spin)
        newSpin2 = c.newSpinEndMove(resn2.pos, newPos2, resn3.pos, resn2.spin)
        
        # Remove residues from lattice
        l.emptyLatticePos(resn1.pos)
        l.emptyLatticePos(resn2.pos)
        
        # Calculate old and new stats
        Moves.oldLocalStats.clean()
        Moves.oldLocalStats.localStats(resn1, l)
        Moves.oldLocalStats.localStats(resn2, l)
        Moves.oldLocalStats.solventStats(newPos1, l)
        Moves.oldLocalStats.solventStats(newPos2, l)
        
        Moves.newLocalStats.clean()
        Moves.newLocalStats.localStats(resn1, newPos1, newSpin1, resn1.fwd, newBkwdn1, resn1.state, resn1.aa, l)
        Moves.newLocalStats.localStats(resn2, newPos2, newSpin2, newFwdn2, resn2.bkwd, resn2.state, resn2.aa, l)
        Moves.newLocalStats.solventStats(resn1.pos, l)
        Moves.newLocalStats.solventStats(resn2.pos, l)
        
        Moves.newLatticeStats.clean()
        Moves.newLatticeStats = l.stats.delta(Moves.newLocalStats, Moves.oldLocalStats)
        
        dE = Moves.newLocalStats.getDeltaE(Moves.oldLocalStats)
        boltz = math.exp(-dE * l.getBetaMoves())
        bb = boltz / (boltz + 1.0)
        
        if l.energyMap is not None:
            l.energyMap.mapStats(l.stats, 1.0 - bb)
            l.energyMap.mapStats(Moves.newLatticeStats, bb)
        
        if dE <= 0 or random.random() < boltz:
            # Accept move
            resn1.spin = newSpin1
            resn2.spin = newSpin2
            resn0.fwd = newFwdn0
            resn1.bkwd = newBkwdn1
            resn2.fwd = newFwdn2
            resn3.bkwd = newBkwdn3
            resn1.pos = newPos1
            resn2.pos = newPos2
            l.setResidue(newPos1, resn1)
            l.setResidue(newPos2, resn2)
            l.stats = Moves.newLatticeStats
            return 1 if dE == 0 else 2
        else:
            # Reject move
            l.setResidue(resn1.pos, resn1)
            l.setResidue(resn2.pos, resn2)
            return 0

    @staticmethod
    def globalMove(l: Lattice, c: Chain) -> int:
        return Moves.translate(l, c) if random.random() < 0.5 else Moves.rotatePoint(l, c)

    @staticmethod
    def translate(l: Lattice, c: Chain) -> int:
        dir = random.randint(0, 5)
        posNew = [Pos() for _ in range(c.N)]
        
        # Check for collisions
        for n in range(c.N):
            res = c.residues[n]
            posNew[n] = res.pos + local[dir]
            posNew[n].periodicBoundary()
            testRes = l.getResidue(posNew[n])
            if testRes is not None and testRes.chainNum != res.chainNum:
                l.sampleCurrentStats()
                return -1
                
        # Calculate old and new stats
        Moves.oldLocalStats.clean()
        Moves.newLocalStats.clean()
        
        # Clear positions from lattice
        for n in range(c.N):
            l.emptyLatticePos(c.residues[n].pos)
            
        for n in range(c.N):
            Moves.oldLocalStats.localStats(c.residues[n], l)
            Moves.oldLocalStats.solventStats(posNew[n], l)
            Moves.newLocalStats.localStats(c.residues[n], posNew[n], l)
            Moves.newLocalStats.solventStats(c.residues[n].pos, l)
            
        Moves.newLatticeStats.clean()
        Moves.newLatticeStats = l.stats.delta(Moves.newLocalStats, Moves.oldLocalStats)
        
        dE = Moves.newLocalStats.getDeltaE(Moves.oldLocalStats)
        boltz = math.exp(-dE * l.getBetaMoves())
        bb = boltz / (boltz + 1.0)
        
        if l.energyMap is not None:
            l.energyMap.mapStats(l.stats, 1.0 - bb)
            l.energyMap.mapStats(Moves.newLatticeStats, bb)
            
        if dE <= 0 or random.random() < boltz:
            # Accept move
            for n in range(c.N):
                c.residues[n].pos = posNew[n]
                l.setResidue(posNew[n], c.residues[n])
            l.stats = Moves.newLatticeStats
            return 1 if dE == 0 else 2
        else:
            # Reject move
            for n in range(c.N):
                l.setResidue(c.residues[n].pos, c.residues[n])
            return 0

    @staticmethod
    def rotatePoint(l: Lattice, c: Chain) -> int:
        rotationalDir = random.randint(0, 5)
        pivot = random.randint(0, c.N - 1)
        
        if c.residues[pivot].state == State.inStrand:
            l.sampleCurrentStats()
            return -1
            
        startRes = 0
        endRes = c.N - 1
        
        if not c.frozen:
            if random.random() < 0.5:
                startRes = pivot
                if pivot < c.N - 1 and c.residues[pivot+1].state == State.inStrand:
                    l.sampleCurrentStats()
                    return -1
            else:
                endRes = pivot
                if pivot > 0 and c.residues[pivot-1].state == State.inStrand:
                    l.sampleCurrentStats()
                    return -1
                    
        if startRes == endRes:
            return -1
            
        # Check for collisions
        posNew = [Pos() for _ in range(c.N)]
        for n in range(startRes, endRes + 1):
            res = c.residues[n]
            posNew[n] = c.rotationPosition(res.pos, rotationalDir, c.residues[pivot].pos)
            testRes = l.getResidue(posNew[n])
            if (testRes is not None and 
                (testRes.chainNum != res.chainNum or 
                 testRes.n < startRes or 
                 testRes.n > endRes)):
                l.sampleCurrentStats()
                return -1
                
        # Check spin directions
        newFwdPivot, newBkwdPivot = Pos(), Pos()
        
        if startRes == pivot:
            newFwdPivot.periodicSubtraction(posNew[pivot+1], posNew[pivot])
            newBkwdPivot = c.residues[pivot].bkwd
            if c.residues[pivot].spin == newFwdPivot:
                l.sampleCurrentStats()
                return -1
        else:
            newBkwdPivot.periodicSubtraction(posNew[pivot-1], posNew[pivot])
            newFwdPivot = c.residues[pivot].fwd
            if c.residues[pivot].spin == newBkwdPivot:
                l.sampleCurrentStats()
                return -1
                
        # Calculate new spins
        newSpins = [Pos() for _ in range(c.N)]
        newFwd = [Pos() for _ in range(c.N)]
        newBkwd = [Pos() for _ in range(c.N)]
        
        for n in range(startRes, endRes + 1):
            newSpins[n] = c.rotationSpin(c.residues[n].spin, rotationalDir)
            newFwd[n] = c.rotationSpin(c.residues[n].fwd, rotationalDir)
            newBkwd[n] = c.rotationSpin(c.residues[n].bkwd, rotationalDir)
            
        newSpins[pivot] = c.residues[pivot].spin
        newFwd[pivot] = newFwdPivot if pivot != c.N - 1 else Pos(0,0,0)
        newBkwd[pivot] = newBkwdPivot if pivot != 0 else Pos(0,0,0)
        
        # Calculate stats
        Moves.oldLocalStats.clean()
        Moves.newLocalStats.clean()
        
        exclStart = startRes + 1 if startRes == pivot else startRes
        exclEnd = endRes - 1 if startRes != pivot else endRes
        
        # Clear lattice positions
        for n in range(exclStart, exclEnd + 1):
            l.emptyLatticePos(c.residues[n].pos)
            
        for n in range(exclStart, exclEnd + 1):
            res = c.residues[n]
            Moves.oldLocalStats.localStats(res, res.pos, l)
            Moves.newLocalStats.localStats(res, posNew[n], newSpins[n], newFwd[n], newBkwd[n], res.state, res.aa, l)
            Moves.oldLocalStats.solventStats(posNew[n], l)
            Moves.newLocalStats.solventStats(res.pos, l)
            
        Moves.newLatticeStats.clean()
        Moves.newLatticeStats = l.stats.delta(Moves.newLocalStats, Moves.oldLocalStats)
        
        dE = Moves.newLocalStats.getDeltaE(Moves.oldLocalStats)
        boltz = math.exp(-dE * l.getBetaMoves())
        bb = boltz / (boltz + 1.0)
        
        if l.energyMap is not None:
            l.energyMap.mapStats(l.stats, 1.0 - bb)
            l.energyMap.mapStats(Moves.newLatticeStats, bb)
            
        if dE <= 0 or random.random() < boltz:
            # Accept move
            for n in range(startRes, endRes + 1):
                c.residues[n].pos = posNew[n]
                l.setResidue(posNew[n], c.residues[n])
                
            for n in range(startRes, endRes + 1):
                c.residues[n].spin = newSpins[n]
                c.residues[n].fwd = newFwd[n]
                c.residues[n].bkwd = newBkwd[n]
                
            l.stats = Moves.newLatticeStats
            return 1 if dE == 0 else 2
        else:
            # Reject move
            for n in range(startRes, endRes + 1):
                l.setResidue(c.residues[n].pos, c.residues[n])
            return 0
