/**
 * \file   Design.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Design
 */


#ifndef _DESIGN_H_
#define _DESIGN_H_

#include <string.h>
#include <stdio.h>
#include <iomanip>
using namespace std;

class Lattice;
class AADistr;

/// Class with methods and objects to design a sequence, given a structure
class Design{
public:
  /// Constructor
  Design();
  /// run a design procedure
  void designProcedure(double beta, long numSteps);
  /// write the statistics of the design procedure to a file
  void writeDesignStats(string file_out);
  /// initialise the Design class, given a PDB structure
  void init(string fn_in);
  void initOldFormat(string fn_in);
  void finalize();

  /// Pointer to lattice containing structure
  Lattice * mainLattice;
  /// Pointer to naturally occuring amino acid distribution
  AADistr * aaDistr;
private:
  void changeAA(int chain);
  void initCountAA();
  double getSystemVar();
  void checkVariance();

};

#endif
