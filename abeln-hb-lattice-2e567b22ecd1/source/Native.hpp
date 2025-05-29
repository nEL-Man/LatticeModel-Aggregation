/**
 * \file   Native.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Native
 */

#ifndef _Native_H_
#define _Native_H_


#include "Lattice.hpp"
#include "Chain.hpp"
#include "Stats.hpp"



/// Class that stores the native contacts of a structure.
/// This is used to calculate the order paremeter native contacts
/// for a given configuration.
class Native{
public:
  /// Constructor, given a Lattice
  Native(Lattice *l);
  /// set native structure - this can consist of multiple chains
  void setNative(Lattice *l);
  /// checks of given pairs of residues is a native contact
  bool isNativeContact(int chainA, int chainB, int resA,  int resB);
  /// check if given energy, is below that of the native structure
  bool isBelowNativeEnergy(int Etot);
  /// gives energy of native contact
  int getNativeEnergy();
  /// checks if given number, is the same as maximum number of native contacts
  bool hasNativeStructure(int Ntot);
private:
  bool ****contacts;
  int totCnat;
  int nativeEnergy;
};



inline
bool Native::isNativeContact(int chainA, int chainB,int resA, int resB){
  return contacts[chainA][chainB][resA][resB];
}

inline
bool Native::isBelowNativeEnergy(int Etot){
  return (Etot < nativeEnergy);
}

inline
bool Native::hasNativeStructure(int Ntot){
  return (Ntot == totCnat);
}


inline
int Native::getNativeEnergy(){
  return nativeEnergy;
}


#endif
