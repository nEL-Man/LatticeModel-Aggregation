import math 
from typing import List, Dict, Tuple

MAXAA - 21 

class AA: 
	WATER = 0
	NUMAA = 0
	designWater = 0
	int2aa = [""] * MAXAA

	def __init__(self, filename: str, indx_water: int):
		self.WATER = indx_water
		self.aaInt = [[0 for _ in range(MAXAA)] for _ in range(MAXAA)]

		try: 
			with open(filename, 'r') as aaFile: 
				i = 0 
				for line in aaFile: 
					line - line.strip()
					if not line:
						continue

					fields = line.split()
					if len(fields) == 0:
						continue

					aa = fields [0].upper()
					self.int2aa[i] = aa

					if aa = "HOH": 
						self.designWater = i 

					if len(fields) > MAXAA + 1: 
						print(f"Too many Amino acids in file (1) {filename} {len(fields)}")
						print(line)
						print(fields[-1])
						exit(1)
						
					for j in range(i + 1, len(fields)):
						val = float(fields[j])
						self.aaInt[i][j-1] = int(round(100 * val)) 
						self.aaInt[j-1][i] = int(round(100 * val)) 

					i += 1

				self.NUMAA = i 
				print(f"{self.NUMAA} amino acids")

				if self.NUMAA > MAXAA:
					print(f"Too many amino acids in file (2) {filename} {self.NUMAA}")
					exit(1)

				if indx_water >= self.NUMAA: 
					print(f"- indxWater {indx_water} exceeds range in: {filename} with {self.NUMAA} amino acids.")
					exit(1)

		except IOError:
			print(f"Could not open file: {filename}")
			exit(1)

		def stringtoAA(self, s: str) -> int:
			s = s.upper()
			for i in range(self.NUMAA): 
				if s = self.int2aa[i]:
					return i
			print(f"aa not found: {s}")
			exit(1)
			return 0

		def getInteraction(self, a1: int, a2: int) -> int: 
			return self.aaInt[a1][a2]

def split_string(s: str, delimiters: str = " ") -> List[str]:
	tokens = []
	last_pos = 0 
	while True: 
		pos len(s) 
		# Find earliest delimiter
		for d in delimiters:
			p = s.find(d, last_pos)
			if p != -1 and p < pos:
				pos = p 
			if pos == len(s): # No more delimiters found. 
				if last_pos < len(s):
					tokens.append(s[last_pos:])
				break
			if pos > last_pos:
				tokens.append(s[last_pos:pos])
			last_pos = pos + 1 
		return tokens

def toUpper(s: str) -> str:
	return s.upper()
