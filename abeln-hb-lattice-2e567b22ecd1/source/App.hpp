/**
 * \file   App.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class App
 */


/*    
    Copyright (C) 2013  Sanne Abeln

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program, see gpl-3.0.txt or see:
    http://www.gnu.org/licenses/.

*/


/*! \mainpage HB-Lattice Documentation
 *
 * \section info_sec Information
 *
 * The header files of the HB-Lattice code have been annotated with
 * Doxygen style comments, on which this documentation is based.
 *  
 * Top level classes are: MonteCarlo (simulation), Design (design) and App (command line options)
 * and the c structured parallel tempering file: ptemp.h and ptemp.cpp
 *
 * \section downl_sec Download
 * 
 * The code is kept in a public git repository on BitBucket.org
 * The latest version of the code may be downloaded from:
 *
 * http://bitbucket.org/abeln/hb-lattice/get/master.tar.gz
 *
 * The repository may be accessed from:
 * 
 * http://bitbucket.org/abeln/hb-lattice
 *
 *
 * \section instal_sec Installation Requirements
 *
 *  To compile this code, Open MPI needs to be installed.
 *  We have compiled the code both on linux or mac like platforms. 
 *  Note that depending on your type of mpi installation, 
 * the finalize command in ptemp.cpp may need to be removed.
 *
 * 
 *
 */








#ifndef _App_Hpp_
#define _App_Hpp_

class MonteCarlo;

#include <string>

using namespace std;

/// Class storing all command line options as class variables
class App{
public:
  static void init(int argc, char * argv[]);
  static MonteCarlo * getMonteCarlo(){return mcSim;};
  static string fn_stats;
  static string fn_eMap ;
  static string fn_hbond ;
  static string fn_in ;
  static string fn_native;
  static string fn_AAdistribution;


  static bool linear_T;
  static bool get_Cext;
  static bool get_PeriodicPDB;
  static bool LAMBDA_COR;
  static bool useOldFormat;
  static bool designProgram;
  static bool AAdistribution;
  static bool checkBelowNative;
  static bool findLowestE;
  static bool writeCluster;
  static bool clustStats;
  static bool nativeStats;
  static bool hbondStats;
  static bool cintStats;
  static bool setAir;
  static bool flipSpins;
  static bool clusterMoves;

  static ofstream moviestream;
  static bool recordMovie;
  static int moviestep;
  static int num_configs;
  static int tableStep;
  static int myRank;

  static string fn_aa;
  static int indx_water;
  
  static double minBeta;
  static double maxBeta;

  static double eBetaMu;

  static double designT;

  // Energy params

  static int HBondE;
  static int StericE;
  static int BBSolvE;

private:
  static MonteCarlo * mcSim;
};




#endif
