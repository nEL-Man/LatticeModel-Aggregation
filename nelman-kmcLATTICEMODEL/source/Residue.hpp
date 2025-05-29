/**
 * \file   Residue.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Residue
 */

#ifndef _Residue_H_
#define _Residue_H_

#include "Pos.hpp"


#define XDIR 0
#define YDIR 1
#define ZDIR 2

/// desribing the state of the Residue
enum State {inCoil,inStrand};

extern Pos local[6];

/// Class containing residue information
class Residue{
public:
  /// amino acid index number
  int aa;
  /// position on the Lattice
  Pos pos;
  /// direction in which the Residue points
  Pos spin;
  /// bkwd direction of the chain, for faster calculations
  Pos bkwd;
  /// fwd direction of the chain, for faster calculations
  Pos fwd;
  /// index of the Residue in the Chain
  int n; 
  /// the index of the Chain
  int chainNum;
  /// the State of the residue, either inStrand or inCoil
  State state;
  /// If masked == true, Cint, Cext, Nint, Next, Hext, Hint statistics are ignored
  bool masked;

  Residue & operator=(const Residue & r);
  Residue(Residue * r);
  Residue(){};
};


inline 
Residue::Residue(Residue * r){
  aa=r->aa;
  pos=r->pos;
  spin = r->spin;
  bkwd=r->bkwd;
  fwd=r->fwd;
  n=r->n;
  chainNum=r->chainNum;
  state=r->state;
  masked = r->masked;
}



inline
Residue & Residue::operator=(const Residue  & r){
  aa=r.aa;
  pos=r.pos;
  spin = r.spin;
  bkwd=r.bkwd;
  fwd=r.fwd;
  n=r.n;
  chainNum=r.chainNum;
  state=r.state;
  masked = r.masked;
  return (*this);
}





#endif
