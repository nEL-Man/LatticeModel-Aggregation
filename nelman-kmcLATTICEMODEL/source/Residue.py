from enum import Enum, auto
from dataclasses import dataclass
from typing import Optional
from pos import Pos  # Assuming you've converted Pos.hpp to pos.py

# Define directions as constants (replacing #define)
XDIR = 0
YDIR = 1
ZDIR = 2

# Enum for residue state (replacing C++ enum)
class State(Enum):
    IN_COIL = auto()
    IN_STRAND = auto()

# Global variable (replacing extern Pos local[6])
local = [Pos() for _ in range(6)]  # Initialize 6 Pos objects

@dataclass
class Residue:
    """Class containing residue information"""
    aa: int                          # amino acid index number
    pos: Pos                         # position on the Lattice
    spin: Pos                        # direction residue points
    bkwd: Pos                        # backward direction
    fwd: Pos                         # forward direction
    n: int                           # index in the Chain
    chainNum: int                    # index of the Chain
    state: State = State.IN_COIL     # default to inCoil
    masked: bool = False             # default to not masked

    # No need for manual copy constructor/assignment operator - 
    # Python handles this automatically with copy.deepcopy() if needed

# Example usage:
if __name__ == "__main__":
    res = Residue(
        aa=1,
        pos=Pos(),
        spin=Pos(),
        bkwd=Pos(),
        fwd=Pos(),
        n=0,
        chainNum=0
    )
    print(res)
