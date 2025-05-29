/**
 * \file   MolBox.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class MolBox
 */


#ifndef _MolBox_hpp_
#define _MolBox_hpp_

#include "Pos.hpp"
#include "Residue.hpp"

#define MAXN 100
#define MAX_CONFIGS 10000
#define EQUIL_STEPS 10000


class Chain;
class Lattice;


/// Struct as helper in MolBox Class to store chain configuration
struct Config{
  Pos positions[MAXN];
  Pos spins[MAXN];
  State states[MAXN];
  int nTotal;
};

/// Class to store and simulate an additional Chain to be inserted via GranCan
class MolBox{
public:
  /// Constructor, given a lattice in which new configurations may be inserted
  MolBox(Lattice * mainLattice);
  /// Sets the chain to be simulated
  void setMolecule(Chain * chain);
  /// create a set of configurations for trial insertions in mainLattice
  void createNewSet(double beta,int num_configs, bool betaChanged);
  void  printConfig(int c);
  /// get the next configuration in set
  Config * getNewConfig();
  /// initialises the interaction matrix
  void setAA(string fn,int water_indx);
  /// get current beta
  double getBeta();
  /// sequence of chain in MolBox
  int sequence[MAXN];
protected:
  /// set of configurations
  Config * configSet[MAX_CONFIGS];
  /// Lattice for creating configurations
  Lattice * box;                       // perhaps need to make this box smaller
  /// chain to be simulated
  Chain * molecule;

private:
  void recordConfig(int insertAt);
  /// equililibration procedure
  void equilibrate(int nsteps);
  /// pointer to currently used configuration
  int currentConfig;
  /// number of Config instances created
  int numConfigs;
 
  /// pointer to mainLattice
  Lattice * mainLattice;
};


#endif






