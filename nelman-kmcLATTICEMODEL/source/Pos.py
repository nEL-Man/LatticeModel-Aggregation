from typing import Tuple, List
import math

LX = 30
LY = 30
LZ = 30

class Pos:
    local = [
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (-1, 0, 0), (0, -1, 0), (0, 0, -1)
    ]
    
    def __init__(self, x: int = 0, y: int = 0, z: int = 0):
        self.x = x
        self.y = y
        self.z = z
    
    def __add__(self, other: 'Pos') -> 'Pos':
        return Pos(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __sub__(self, other: 'Pos') -> 'Pos':
        return Pos(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def __neg__(self) -> 'Pos':
        return Pos(-self.x, -self.y, -self.z)
    
    def __eq__(self, other: 'Pos') -> bool:
        return self.x == other.x and self.y == other.y and self.z == other.z
    
    def __ne__(self, other: 'Pos') -> bool:
        return not (self == other)
    
    def __getitem__(self, index: int) -> int:
        return [self.x, self.y, self.z][index]
    
    def __setitem__(self, index: int, value: int):
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        elif index == 2:
            self.z = value
    
    def periodicBoundary(self):
        self.x = self.x % LX
        self.y = self.y % LY
        self.z = self.z % LZ
    
    def periodicSubtraction(self, p1: 'Pos', p2: 'Pos'):
        self.x = p1.x - p2.x
        self.y = p1.y - p2.y
        self.z = p1.z - p2.z
        
        if self.x > (LX // 2):
            self.x -= LX
        elif self.x < -(LX // 2):
            self.x += LX
        
        if self.y > (LY // 2):
            self.y -= LY
        elif self.y < -(LY // 2):
            self.y += LY
        
        if self.z > (LZ // 2):
            self.z -= LZ
        elif self.z < -(LZ // 2):
            self.z += LZ
    
    def printCout(self):
        print(f"x:{self.x} y:{self.y} z:{self.z}")
    
    def toString(self) -> str:
        return f"x:{self.x} y:{self.y} z:{self.z}"
    
    @staticmethod
    def orthogonal(p1: 'Pos', p2: 'Pos') -> bool:
        return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z) == 0

# Non-member functions
def getAngleUnitV(p1: Pos, p2: Pos) -> int:
    dot = p1.x * p2.x + p1.y * p2.y + p1.z * p2.z
    return math.acos(dot)

def getAnglePositionV(p1: Pos, p2: Pos, p3: Pos) -> int:
    v1 = p2 - p1
    v2 = p3 - p2
    return getAngleUnitV(v1, v2)

def pickNeighbour90(p1: Pos, p3: Pos, random_n: int) -> Pos:
    # Implementation depends on specific requirements
    pass

def pickNeighbour120(p1: Pos, p2: Pos, p3: Pos) -> Pos:
    # Implementation depends on specific requirements
    pass
