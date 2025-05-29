/**
 * \file   Chain.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Chain
 */


#ifndef _Chain_H_
#define _Chain_H_

#include "Residue.hpp"
#include "Stats.hpp"


/// defines the maximum number of residues held by the chain
#define MAX_RES 100


class Lattice;

/// Class containing the information about a peptide chain
/// and contains Residue instances
class Chain{
public:
  /// constructor, given a lattice, the chain needs to be inserted on, and id of the chain
  Chain(Lattice *l,int chnum);
  /// destructor: deletes the residues
  ~Chain();

  void setChainNum(int n);
  void setFwdBkwd(Residue * res);
  void setFwdBkwdSpins();
  bool checkFwdBkwd(Residue * res);
  bool checkAllFwdBkwd();
  
  int getCext();
  bool hasCext();
  void setRandomSpins();
  void newRandomSpin(Residue * res);

  void setStateResidues();

  //// helper functions for rotating chains
  int rotationCoord(int x,int y,int z, int dir, int rotationDir,int pivot);
  Pos rotationSpin(Pos spinOld, int rotationDir);
  Pos rotationPosition(Pos posOldi, int rotationDir,Pos posPivot);
  bool strandPossible(Residue * res);

  //// helper functions for spins
  Pos newSpinEndMove(Pos oldPos, Pos newPos, Pos posTail, Pos oldSpin);
  Pos newSpinPosCornerFlip(Pos oldSpin,Pos n1, Pos n2); 


  /// size of the chain
  int N;
  /// array containing pointers to Residues in the Chain
  Residue * residues[MAX_RES];
  bool frozen;
  bool locked;

  int chainNum;

private:
  Lattice * l;
 

  // Variables to be used in calculation of new stats
  Stats oldLocalStats;
  Stats newLocalStats;
  Stats newLatticeStats;

 
};


#endif
