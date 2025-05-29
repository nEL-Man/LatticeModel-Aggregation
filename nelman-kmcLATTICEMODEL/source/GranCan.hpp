/**
 * \file   GranCan.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class GranCan
 */


#ifndef _GranCan_H_
#define _GranCan_H_

#include <string.h>
#include <iostream>


using namespace std;


/// if defined does not simulate empty box explicityly
#define CORRECTION_0CHAINS

class MolBox;
class Lattice;

/// Class to enable Grand Canonical simulations.
/// It is part of the MonteCarlo simulation
/// it contains an instance of MolBox
/// Its methods are used in Moves and MonteCarlo
class GranCan{
public:
  GranCan(Lattice * mainLattice);
 
  int trialInsertChain();
  int trialDeleteChain();
  int getNumFreeChains();
  int clusterMove();
  int getCext();

  /// simulation box for generating configurations for trial insertions
  MolBox * molbox;
  double pGlobal;
  double pInsDel;
  double pChangeStrandCoil;
  double pClusterMove;

  int steps;
  int total_swaps;
  double totalEmptySteps;
  double Volume;

private:
  Lattice * mainLattice;


};


#endif
