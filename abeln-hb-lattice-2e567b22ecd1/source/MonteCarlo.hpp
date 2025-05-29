/**
 * \file   MonteCarlo.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class MonteCarlo
 */


#ifndef _MonteCarlo_H_
#define _MonteCarlo_H_

#include <string.h>
#include <iostream>
using namespace std;


class MolBox;
class Lattice;
class GranCan;

/// Class containing methods and object to run the Monte Carlo simulation.
/// It contains instances of Lattice and GranCan.
/// It is called by ptemp
class MonteCarlo{
public:
  /// Constructor
  MonteCarlo();
  /// Run the MonteCarlo simulation
  int MC(double beta, int betaId, bool betaChanged);
  /// Get the statistics of a MC run
  void getMCStats(int num_procs, int myRank);
  /// get external contacts
  int getCext();

  /// Lattice on which the simulation is performed
  Lattice * mainLattice;
  /// GranCan for grand canonical simulations
  GranCan * grancan;
  double pGlobal;
  double pInsDel;
  double pChangeStrandCoil;
  double pClusterMove;
  /// number of steps for each MC run
  int steps;
  /// number of swaps in parallel tempering, see ptemp.cpp
  int total_swaps;


};





#endif
