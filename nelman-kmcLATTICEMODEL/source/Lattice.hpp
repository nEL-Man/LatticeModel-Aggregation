/**
 * \file   Lattice.hpp
 * \date   May 2013
 * \author SA
 * \brief  Defines class Lattice
 */


#ifndef _Lattice_H_
#define _Lattice_H_

/// the maximum number of chains on a lattice
#define MAX_CHAINS 50000


#include "Residue.hpp"
#include "Stats.hpp"
#include "EnergyMap.hpp"

#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string.h>

using namespace std;

class Chain;
class Stats;
class Native;
class EnergyMap;
class AA;
class MolBox;
struct Config;

enum Solvent {water,air};

/// Class containing the lattice cells on which a Chain is simulated
class Lattice{
public:
  /// Constructor
  Lattice();

  Residue * getResidue(const Pos & p);
  void emptyLatticePos(const Pos & p);
  void setResidue(const Pos &  p,Residue * r);
  
  Solvent getSolvent(const Pos & p);

  void setSolvent();
  void readPDBMultiChain(string s);
  void oldReadPDBMultiChain(string s);
  /// Writes current configuration on lattice to given filename
  void writePDBMultiChain(string s);
  /// appends current configuration PDB movie file
  void writePDBMultiChain(ofstream & outs,int moviestep);
  void printPeriodicPDB();
  int insertChain(Config * config, int * sequence, double beta);
  int deleteChain(int chainNum);
  void setAA(string fn, int water_indx);
  bool checkStats();
  void resetStats();
  double getBetaMoves();
  double getRealBeta();
  void setBetaMoves(double b);
  /// Sets native structure, given filename of a PDB structure in lattice format
  void setNative(string fn_native);
  void applyPeriodicBoundaries();
  void sampleCurrentStats();
  bool checkForClashesAndTouches(const Pos & p);

  /// array containing all chains on the Latticee
  Chain * chains[MAX_CHAINS];
  /// the lattice 
  Residue * r[LX][LY][LZ];
  /// contains global statistics about the lattice
  Stats stats;
  /// the number of chains on the Lattice
  int nChains;
  /// pointer to a native structure
  Native * native;
  /// pointer to EnergyMap: mapping energy statistics of the simulation
  EnergyMap * energyMap;
  /// pointer the amino acid interaction matrix
  AA * aaInt;
  /// dummy pointer to indicate AIR interface
  Residue * AIR;
private:
  double betaMoves;

  Solvent solv[LX][LY][LZ];

  //int freeChains;
 
 
};




inline
Residue * Lattice::getResidue(const Pos & p){
  return r[p.x][p.y][p.z];
}


inline
void Lattice::emptyLatticePos(const Pos & p){
  r[p.x][p.y][p.z]=NULL;
}

inline
void Lattice::setResidue(const Pos & p,Residue * res){
  r[p.x][p.y][p.z]=res;
}


inline
Solvent Lattice::getSolvent(const Pos & p){
  return solv[p.x][p.y][p.z];
}

inline
void Lattice::setBetaMoves(double m){
  betaMoves=m;
}


inline
double Lattice::getBetaMoves(){
  return betaMoves;
}

inline
double Lattice::getRealBeta(){
  return 100*betaMoves;
}

inline
void Lattice::sampleCurrentStats(){
  if(energyMap){
    energyMap->mapStats(stats,1.0);
  }
}


inline 
bool Lattice::checkForClashesAndTouches(const Pos & p){
  if( r[p.x][p.y][p.z] !=NULL){
    return true;
  }
  for(int k=0;k<6;k++){
	Pos posNB = local[k] + p;
	posNB.periodicBoundary();
	if(getResidue(posNB)!=NULL){
	  return true;
	}
  }
  return false;
}


#endif
