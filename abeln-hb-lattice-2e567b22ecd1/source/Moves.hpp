/**
 * \file   Moves.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Moves
 */


#ifndef _Moves_H_
#define _Moves_H_

#include "Residue.hpp"
#include "Stats.hpp"


class Lattice;

/// Class containing methods to generate moves in the MonteCarlo simulation
class Moves{
public:
  static int localMove(Lattice *l, Chain * c);
  static int globalMove(Lattice *l,Chain * c);
  static int rotatePoint(Lattice *l,Chain * c);
  static int changeStrandCoil(Lattice *l,Chain * c);
  static int shuffleSpin(Lattice *l,Chain *c);
  static int crankshaft(Lattice *l,Chain *c,int n1, int n2,int px, int py, int pz);
  static int translate(Lattice *l,Chain *c);

private:



  // Variables to be used in calculation of new stats
  static Stats oldLocalStats;
  static Stats newLocalStats;
  static Stats newLatticeStats;
  //  double betaMoves;
 
};


#endif
