/**
 * \file   AA.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class AA
 */



#ifndef _AA_H_
#define _AA_H_

#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <iostream>

#define MAXAA 21

using namespace std;

/// Class containing the interaction matrix for Amino Acids
class AA{
public:
  /// constructor, given a filename and index for the solvent row
  AA(string filename,int water_indx);
  /// Gives the interaction between a1 and a2
  int getInteraction(int a1,int a2);
  /// 3 letter amino acid code to single capital code
  int stringtoAA(string s);
  static string int2aa[MAXAA];
  /// total number of amino acids in matrix 
  static int NUMAA;
  /// index of water
  static int WATER;
  /// index for water in design procedure 
  static int designWater;
private:
  /// matrix containing interactions
  int aaInt[MAXAA][MAXAA];
};


inline
int AA::getInteraction(int a1,int a2){
  return aaInt[a1][a2];
}

#endif
