from typing import Optional, List
import argparse
import sys
from MonteCarlo import MonteCarlo

class App: 
  fn_aa = "aa.txt"
  indx_water = -1
  minBeta = 0.01
  maxBeta = 0.05 
  fn_stats + "stats_" 
  fn_eMap = "eMap_" 
  fn_hbond = "hbond_"
  fn_in = "RWWLY10_solo.pdb"
  fn_native = ""
  fn_AAdistribution = ""

  # Flags
  nativeStats = False
  linear_T = False
  get_Cext = False
  get_PeriodicPDB = False 
  LAMBDA_COR = False
  useOldFormat = False
  designProgram = False
  AAdistribution = False
  checkBelowNative = False
  findLowestE = False 
  writeCluster = False
  hbondStats = False
  cintStats = False
  clustStats = False
  setAir = False
  flipSpins - True
  clusterMoves = False

  eBetaMu = -0.0000005
  recordMovie = False
  moviestep = 0 
  num_configs = 10000
  tableStep = 0 
  myRank = -1

  # Energy Terms
  HBondE = -50
  StericE = 55
  BBSolveE = 0

  designT = 0.0
  mcSim: Optional[MonteCarlo] = None
  
@classmethod
def init(cls, argv: List[str]): 
  cls.mcSim = MonteCarlo()
  minT = 0.0
  maxT = 0.0

  parser = argparse.ArgumentParser(description='KMC-Lattice Simulation Parameters')
  parser.add_argument('-minB', type=float, help='Minimum beta Value')
  parser.add_argument('-maxB', type=float, help='Maximum beta value')
  parser.add_argument('-minT', type=float, help='Minimum Temperature')
  parser.add_argument('-maxT', type=float, help='Maximum Temperature')
  parser.add_argument('-S', type=int, help='Simulation steps')
  parser.add_argument('-W', type=int, help='Total swaps')
  parser.add_argument('-i', type=str, help='Input file')
  parser.add_argument('-oldFormat', action='store_true', help='Use old format')
  parser.add_argument('-ebmu', type=float, help='eBetaMu value')
  parser.add_argument('-aa', type=str, help='Amino acid file')
  parser.add_argument('-indxWater', type=int, help='Water index')
  parser.add_argument('-setAir', action='store_true', help='Set air')
  parser.add_argument('-native', type=str, help='Native structure file')
  parser.add_argument('-lT', action='store_true', help='Linear temperature')
  parser.add_argument('-noSwap', action='store_true', help='Disable swapping')
  parser.add_argument('-getCext', action='store_true', help='Get Cext')
  parser.add_argument('-getPeriodicPDB', action='store_true', help='Get periodic PDB')
  parser.add_argument('-lambdaCor', action='store_true', help='Lambda correction')
  parser.add_argument('-movie', action='store_true', help='Record movie')
  parser.add_argument('-checkBelowNative', action='store_true', help='Check below native')
  parser.add_argument('-findLowestE', action='store_true', help='Find lowest energy')
  parser.add_argument('-writeCluster', action='store_true', help='Write cluster')
  parser.add_argument('-clustStats', action='store_true', help='Cluster statistics')
  parser.add_argument('-hbondStats', action='store_true', help='HBond statistics')
  parser.add_argument('-CintStats', action='store_true', help='Cint statistics')
  parser.add_argument('-AAdistr', type=str, help='AA distribution file')
  parser.add_argument('-noSpinFlips', action='store_true', help='Disable spin flips')
  parser.add_argument('-clusterMoves', action='store_true', help='Cluster moves')
  parser.add_argument('-designT', type=float, help='Design temperature')
  parser.add_argument('-HBondE', type=int, help='HBond energy')
  parser.add_argument('-StericE', type=int, help='Steric energy')
  parser.add_argument('-BBSolvE', type=int, help='BBSolv energy')

args = parser.parse_args(argv[1:])
        
  if args.minB:
    cls.minBeta = args.minB
  if args.maxB:
    cls.maxBeta = args.maxB
  if args.minT:
    minT = args.minT
  if args.maxT:
    maxT = args.maxT
  if args.S:
    cls.mcSim.steps = args.S
  if args.W:
    cls.mcSim.total_swaps = args.W
  if args.i:
    cls.fn_in = args.i
  if args.oldFormat:
    cls.useOldFormat = True
  if args.ebmu:
    cls.eBetaMu = args.ebmu
  if args.aa:
    cls.fn_aa = args.aa
  if args.indxWater:
    cls.indx_water = args.indxWater
  if args.setAir:
    cls.setAir = True
  if args.native:
    cls.fn_native = args.native
    cls.nativeStats = True
  if args.lT:
    cls.linear_T = True
  if args.noSwap:
    cls.swapT = False
  if args.getCext:
    cls.get_Cext = True
  if args.getPeriodicPDB:
    cls.get_PeriodicPDB = True
  if args.lambdaCor:
    cls.LAMBDA_COR = True
  if args.movie:
    cls.recordMovie = True
  if args.checkBelowNative:
    cls.checkBelowNative = True
  if args.findLowestE:
    cls.checkBelowNative = True
    cls.findLowestE = True
  if args.writeCluster:
    cls.writeCluster = True
  if args.clustStats:
    cls.clustStats = True
  if args.hbondStats:
    cls.hbondStats = True
  if args.CintStats:
    cls.cintStats = True
  if args.AAdistr:
    cls.AAdistribution = True
    cls.fn_AAdistribution = args.AAdistr
  if args.noSpinFlips:
    cls.flipSpins = False
  if args.clusterMoves:
    cls.clusterMoves = True
  if args.designT:
    cls.designProgram = True
    cls.designT = args.designT
  if args.HBondE:
    cls.HBondE = args.HBondE
  if args.StericE:
    cls.StericE = args.StericE
  if args.BBSolvE:
    cls.BBSolvE = args.BBSolvE
      
  if minT > 0:
    cls.maxBeta = 0.01 / minT
    print(f"minBeta {cls.maxBeta}")
  if maxT > 0:
    cls.minBeta = 0.01 / maxT
    print(f"maxBeta {cls.minBeta}")
      
  if cls.designProgram:
    cls.swapT = False
      # Design implementation would go here
      
  cls.mcSim.mainLattice.energyMap.initTableStats()
  
  cls.num_configs = int(cls.mcSim.steps * cls.mcSim.pInsDel * 1.1)
  if cls.num_configs > cls.MAX_CONFIGS:
    cls.num_configs = cls.MAX_CONFIGS
    print(f"num_configs per set: {cls.num_configs}")
  
  if cls.eBetaMu < 0:
    cls.mcSim.pInsDel = 0.0
      
  cls.mcSim.grancan.molbox.setAA(cls.fn_aa, cls.indx_water)
  cls.mcSim.mainLattice.setAA(cls.fn_aa, cls.indx_water)
  
  if not cls.useOldFormat:
    cls.mcSim.mainLattice.readPDBMultiChain(cls.fn_in)
  else:
    cls.mcSim.mainLattice.oldReadPDBMultiChain(cls.fn_in)
      
  if cls.nativeStats:
    if cls.fn_native == "":
      cls.fn_native = cls.fn_in
      cls.mcSim.mainLattice.setNative(cls.fn_native)
      print("set native")
      
  print(f"read file {cls.fn_in}")
  
  if cls.setAir:
    cls.mcSim.mainLattice.setSolvent()
    cls.mcSim.mainLattice.resetStats()
      
  if cls.eBetaMu > 0:
    cls.mcSim.grancan.molbox.setMolecule(cls.mcSim.mainLattice.chains[0])
    print("molbox molecule set")
      
  if cls.get_Cext:
    cls.mcSim.getCext()
    sys.exit(0)
  if cls.get_PeriodicPDB:
    cls.mcSim.mainLattice.printPeriodicPDB()
    sys.exit(0)
      
  if cls.recordMovie:
    cls.moviestream = open("movie.pdb", 'w')
      
  # Initialize and print stats
  newStats = Stats()
  newStats.getLatticeStats(cls.mcSim.mainLattice)
  print("initial stats:")
  newStats.printCout()
