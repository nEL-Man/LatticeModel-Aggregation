/**
 * \file   Stats.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Stats
 */


#ifndef _Stats_H_
#define _Stats_H_

#include "Residue.hpp"

//#include "EnergyMap.hpp"



#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string.h>


//#define StrandE 0
//#define StericE 55
//#define HBondE -50

//#define StericE 100
//#define HBondE -80

//#define CROSSNBS
//#define CROSS_FACTOR 2

//#define EXCL_THR_C

using namespace std;

class Lattice;
class Chain;

/// Class containing simulation statistics.
/// This is used to get changes in order parameters before and after a move.
/// Stats is used by Lattice, Moves, MonteCarlo, Design
/// Many of its methods are defined inline in StatsMethods.hpp

class Stats{
public:
  //inline StatsMethods.hpp
  Stats(); 
  Stats operator+(const Stats & s) const;
  Stats & operator=(const Stats & s);
  Stats & operator+= (const Stats& s);
  Stats & operator-= (const Stats& s);
  const bool  operator!= (const Stats& s) const;
  void clean();
  int getCext(){return Cext;};
  int getEtot(){return Etot;};

  Stats  delta(const Stats & Snew, const Stats & Sold)const;
  int getDeltaE(const Stats & old)const;  
  
  void localStats(Residue * res, Lattice *l);
  void localStats(Residue * res, Pos p, Lattice *l);
  void localStats(Residue * res,const Pos & pos, const Pos & spin,State state,int resAA, Lattice * l);
  void localStats(Residue * res,const Pos & pos, const Pos & spin, const Pos & fwd, const Pos & bkwd, State state,int resAA, Lattice * l);
  void localStats(Residue * res,const Pos & pos,const  Pos & spin, State state, Lattice * l);

  void localStatsExclude(Residue * res,const Pos & pos,const  Pos & spin,State state,  Lattice * l,int start,int end);
  
  void solventStats(Pos & pos,Lattice *l);
  //Stats.cpp
  void getLatticeStats(Lattice * l);
  void get_Eint_Cint(Chain*  c,Lattice *l);
  void old_get_Eint_Cint(Chain*  c, Lattice *l);
  void printCout();
  string print2string();
  
  //static void setLattice(Lattice * l);
#ifdef CROSSNBS
  static void setLoopUpCrossNbs();
  static Pos* getCrossNbsSpinDir(Pos spin);
#endif
private:
  int Etot;
  int Ctot;  
  int Cext;
  int Eint;
  int Cint;
  int Eext;
  int Esol;
  friend class EnergyMap;
  
  int Nint;
  int Next;
  int Ntot;

  int Hint;
  int Hext;
  int Htot;

#ifdef CROSSNBS
  static Pos loopUpCrossNbs[3][3][3][4];
#endif
  //static Lattice * l;
};

/*** NON-MEMBER FUNCTIONS ***/





#endif
