/**
 * \file   AADistr.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class AADistr
 */


#ifndef _AADISTR_HPP_
#define _AADISTR_HPP_

#include <string.h>
#include <stdio.h>
#include <iomanip>

#include "AA.hpp"

using namespace std;

class Lattice;


/// Class containing natural occurring amino acid distribution
/// to aid with the Design process
class AADistr{
public:
  /// constructor
  AADistr();
  /// initialisation of AADistr
  /// @param[in] filename containing the distribution
  /// @param[in] aaInt - the already initialised interaction matrix
  void   init(string filename, AA * aaInt);
  
  /// check composition of a sequence
  void setSequence(Lattice *l);
  /// get distance between sequence and natural occurring distribution
  double getDistance();
  /// get distance after subsitution of amino acids
  double getDistanceSubst(int aa1, int aa2);
  
  /// updates freqency table
  void substitute(int aa1, int aa2,double dist);

  /// send current distribution to string
  string toString();
  bool check();

private:
  /// frequencies as observed in PDB (from file)
  double freqPDB[MAXAA];
  /// expected frequencies given the sequence length
  double freqExp[MAXAA];
  /// frequencies in current sequence
  double freqSeq[MAXAA];
  void setLength(int length);
  int seqLength;
  int numaa;
  AA * aaInt;
  double systemDistance;
  double sumFreq;
};


inline
void AADistr::substitute(int aa1, int aa2, double dist){
  freqSeq[aa1]--;
  freqSeq[aa2]++;
  systemDistance += dist;
}

inline
double AADistr::getDistanceSubst(int aa1, int aa2){
  if(aa1==aa2)return 0;
  return 2*((freqExp[aa1] - freqSeq[aa1])- (freqExp[aa2] - freqSeq[aa2]-1.0)); 
}

#endif
