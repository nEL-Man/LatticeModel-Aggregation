from typing import List 
from AA import AA, MAXAA 

class AADistr:
	def__init__(self):
		self.freqPDB = [0.0] * MAXAA
		self.freqExp = [0.0] * MAXAA
		self.freqSeq = [0.0] * MAXAA 
		self.seqLength = 0 
		self.numaa = 0 
		self.aaInt = None
		self.systemDistance = 0.0
		self.sumFreq = 0 

	def init(self, filename: str, aai: AA): 
		self.aaInt = aai
		self.numaa = AA.NUMAA 
		self.sumFreq = 0

		try:
			with open(filename, 'r') as aaFile: 
				i = 0
				for line in aaFile: 
					line = line.strip()
					if not line: 
						continue

					parts = line.split()
					if len(parts) < 2: 
						continue

					aa = parts[0].upper()
					freq = int(parts[1])
					print(aa)

					if aa != "HOH":
						if self.aaInt.stringtoAA(aa) != i:
							print(f"Amino acid order in file {filename} does not match interaction matrix")
							exit(1)
						self.freqPDB[i] = freq
						self.sumFreq += freq
					i += 1
		except IOError: 
			print(f"Could not open file: {filename}")
			exit(1)

	def setLength(self, length: int):
		self.seqLength = length
		for i in range(self.numaa):
			if i != AA.WATER:
				self.freqExp[i] = self.freqPDB[i] * self.seqLength / self.sumFreq
			else: self.freqExp[i] = 0 

	def setSequence(self, 1): 
		length = 0 
		for i in range(AA.NUMAA):
			self.freqSeq[i] = 0

		for nc in range(l.nChains): 
			ch = l.chains[nc]
			for n in range(ch.N): 
				res = ch.residues[n]
				self.freqSeq[res.aa] += 1
				length += 1

		self.setLength(length)
		self.systemDistance = self.getDistance()

	def getDistance(self) -> float:
		dist = 0.0 
		for i in range(self.numaa):
			if i != AA.WATER:
				dist += self.freqSeq[i] - self.freqExp[i]) ** 2 

		return dist

	def check(self) -> bool:
		calcD = self.getDistance()
		if calcD = self.systemDistance:
			print(f"ERROR in design process accumulative distance: {self.systemDistance} does not match calc distance: {calcD}")
			return False
		return True

	def substitute(self, aa1: int, aa2: int, dist: float):
		self.freqSeq[aa1] -= 1	
		self.freqSeq[aa2] += 1 
		self.systemDistance += dist

	def getDistanceSubst(self, aa1: int, aa2: int) -> float: 
		if aa1 = aa2:
			return 0
		return 2 * ((self.freqExp[aa1] - self.freqSeq[aa1]) - (self.freqExp[aa2] - self.freqSeq[aa2] - 1.0))

	def toString(self -> str: 
		import io 
		output = io.StringIO()

		output.write(f"System distance {self.systemDistance}\n") 
		output.write("SEQUENCE\n")
		1 = 0
		for aa in range(AA.NUMAA): 
			output.write(f"{AA.int2aa[aa]} {self.freqSeq[aa]}, ")
			1 += self.freqSeq[aa]
		output.write(f" Length {l}\n")

		output.write("PDB\n")
		checksum = 0 
		for aa in range(AA.NUMAA): 
			rounded = round(self.freqExp[aa], 1)
			output.write(f"{AA.int2aa[aa]} {rounded}, ")
			checksum += self.freqExp[aa]
		output.write(f"SUM: {checksum}\n"

		return output.getValue()
			     
