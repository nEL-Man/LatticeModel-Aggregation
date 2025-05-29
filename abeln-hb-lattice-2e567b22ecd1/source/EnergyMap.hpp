/**
 * \file   EnergyMap.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class EnergyMap
 */


#ifndef _EnergyMap_H_
#define _EnergyMap_H_

#include <string.h>
//#include <ext/hash_map>
#include <tr1/unordered_map>
#include <fstream>

#include "ptemp.h"
#include "Stats.hpp"

#define QMAX 10000

/// struct containing information on largest cluster
struct ClusterInfo{
  int clusters;
  int largestCluster;
  int totalChains;
};



using namespace std;

class Stats;
class Lattice;

/// Class containing a hash of sampled energies during simulation.
/// It is stored on the Lattice and
/// called from the MonteCarlo simulation
class EnergyMap{
public:
  EnergyMap(Lattice * l);
  void setBetaID(int bID);
  void printGnuPlotDataWR(string fn_stats,int num_procs, int myRank);
  void printEMap(string fn_stats,int num_procs, int myRank);
  void mapStats(Stats &s, double boltz);
  void printHBondMap(string fn_hbond,int num_procs, int myRank);
  void initTableStats();
  void writeTableStats(Stats & stats, int step,double weight, ClusterInfo ci);
  void finalizeTableStats();

private:
  void checkBelowNative(Stats &s, double boltz);
  double energyMap[MAX_PROCS][QMAX];
  double hbondMap[MAX_PROCS][QMAX];
  // __gnu_cxx::hash_map<int,double> eMap[MAX_PROCS];
  tr1::unordered_map<int,double> eMap[MAX_PROCS];
  int betaID;
  Lattice *l;
  ofstream tableStatsStream;
};




#endif
