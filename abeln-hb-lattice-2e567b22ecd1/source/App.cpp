#include "App.hpp"
#include "GranCan.hpp"
#include "MonteCarlo.hpp"
#include "Design.hpp"
#include "ptemp.h"
#include "MolBox.hpp"
#include "Lattice.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"

#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>



using namespace std;

string App::fn_aa ="aa.txt";
int App::indx_water=-1;

double App::minBeta=0.01;
double App::maxBeta=0.05;

string App::fn_stats="stats_";
string App::fn_eMap ="eMap_";
string App::fn_hbond ="hbond_";

string App::fn_in = "RWWLY10_solo.pdb";
string App::fn_native="";
string App::fn_AAdistribution="";


bool App::nativeStats=false;
bool App::linear_T=false;
bool App::get_Cext=false;
bool App::get_PeriodicPDB=false;
bool App::LAMBDA_COR=false;
bool App::useOldFormat=false;
bool App::designProgram=false;
bool App::AAdistribution=false;
bool App::checkBelowNative=false;
bool App::findLowestE=false;
bool App::writeCluster=false;
bool App::hbondStats=false;
bool App::cintStats=false;
bool App::clustStats=false;

bool App::setAir=false;

bool App::flipSpins=true;
bool App::clusterMoves=false;

double  App::eBetaMu = -0.0000005;

ofstream App::moviestream;
bool App::recordMovie=false;
int App::moviestep=0;
int App::num_configs= 10000;
int App::tableStep=0;

int App::myRank=-1;


/// Energy terms
int App::HBondE=-50;
int App::StericE=55;
int App::BBSolvE=0;


double App::designT;

MonteCarlo * App::mcSim=NULL;



void App::init(int argc,char *argv[]){
  mcSim= new MonteCarlo();
  
  double minT=0,maxT=0;

  int optind=1;
#ifdef CROSSNBS
  Stats::setLoopUpCrossNbs();
#endif
  while ((optind < argc) && (argv[optind][0]=='-')) {
    string sw = argv[optind];
    //cout<<"* "<< sw <<endl;
    if (sw=="-minB") {
      optind++;
      minBeta = atof(argv[optind]);
    }else if(sw=="-maxB") {
      optind++;
      maxBeta = atof(argv[optind]);
    }else if (sw=="-minT") {
      optind++;
      minT = atof(argv[optind]);
    }else if(sw=="-maxT") {
      optind++;
      maxT = atof(argv[optind]);
    }else if(sw=="-S"){
      optind++;
      mcSim->steps = atoi(argv[optind]);
    }else if(sw=="-W"){
      optind++;
      mcSim->total_swaps = atoi(argv[optind]);
    }else if(sw=="-i"){
      optind++;
      fn_in.assign(argv[optind]);
    }else if(sw=="-oldFormat"){
      optind++;
      useOldFormat=true;
      fn_in.assign(argv[optind]);
    }else if (sw=="-ebmu") {
      optind++;
      eBetaMu = atof(argv[optind]);
    }else if(sw=="-aa"){
      optind++;
      cout<< argv[optind]<<endl;
      fn_aa.assign(argv[optind]);
    }else if(sw=="-indxWater"){
      optind++;
      indx_water= atoi(argv[optind]);
    }else if(sw=="-setAir"){
      setAir=true;
    }else if(sw=="-native"){
      optind++;
      fn_native.assign(argv[optind]);
      nativeStats=true;
    }else if(sw=="-lT"){
      linear_T=true;
    }else if(sw=="-noSwap"){
      swapT=false;
    }else if(sw=="-getCext"){
      get_Cext=true;
    }else if(sw=="-getPeriodicPDB"){
      get_PeriodicPDB=true;
    }else if(sw=="-lambdaCor"){
      LAMBDA_COR=true;
    }else if(sw=="-movie"){
      recordMovie=true;
    }else if(sw=="-checkBelowNative"){
      checkBelowNative=true;
    }else if(sw=="-findLowestE"){
      checkBelowNative=true;
      findLowestE=true;
    }else if(sw=="-writeCluster"){
      writeCluster=true;
    }else if(sw=="-clustStats"){
      clustStats=true;
    }else if(sw=="-hbondStats"){
      hbondStats=true;
    }else if(sw=="-CintStats"){
      cintStats=true;
    }else if(sw=="-AAdistr"){
      AAdistribution=true;
      optind++;
      fn_AAdistribution.assign(argv[optind]);
    }else if(sw=="-noSpinFlips"){
      flipSpins=false;
    }else if(sw=="-clusterMoves"){
      clusterMoves=true;
    }else if(sw=="-designT"){
      designProgram=true;
      optind++;
      designT = atof(argv[optind]);
    }else if(sw=="-HBondE"){
      optind++;
      HBondE = atoi(argv[optind]);
    }else if(sw=="-StericE"){
      optind++;
      StericE = atoi(argv[optind]);
    }else if(sw=="-BBSolvE"){
      optind++;
      BBSolvE = atoi(argv[optind]);
    }else{
      cout << "unknown option:"<< sw<<":"<<endl;
      exit(1);
    }
    optind++;
  } 
 
  if(minT >0){maxBeta= 0.01/minT;cout<<"minBeta "<<maxBeta<<endl;}
  if(maxT >0){minBeta= 0.01/maxT;cout<<"maxBeta "<<minBeta<<endl;}
 
  if(designProgram){
    swapT=false;
    //Design * des= new Design();
    //if(! useOldFormat){
    //  des->init(fn_in);
    //}else{
    //  des->initOldFormat(fn_in);
    //}
    //cout<<"b"<<endl;
    //des->designProcedure((minBeta + maxBeta)/2.0, mcSim->steps);
    //exit(0);
  }

  mcSim->mainLattice->energyMap->initTableStats();


  num_configs= (int)(mcSim->steps * mcSim->pInsDel * 1.1);
  if(num_configs > MAX_CONFIGS)num_configs =MAX_CONFIGS;
  cout<< "num_configs per set: "<<  num_configs<< endl;
  if(eBetaMu<0){mcSim->pInsDel=0.0;}
 
  mcSim->grancan->molbox->setAA(fn_aa,indx_water);
  mcSim->mainLattice->setAA(fn_aa,indx_water);
  if(! useOldFormat){
    mcSim->mainLattice->readPDBMultiChain(fn_in);
  }else{
    mcSim->mainLattice->oldReadPDBMultiChain(fn_in);
  }

  if(nativeStats){
     if(fn_native== ""){
       fn_native = fn_in;
     }
     mcSim->mainLattice->setNative(fn_native);
     cout<< "set native "<<endl;
  }

  cout<< "read file "<<fn_in<<endl;
  
  if(setAir){
    mcSim->mainLattice->setSolvent();
    //recalculate stats
    mcSim->mainLattice->resetStats();
  }

  if(eBetaMu>0){
    mcSim->grancan->molbox->setMolecule(mcSim->mainLattice->chains[0]);
    cout<< "molbox molecule set "<<endl;
  }
 
  
  if(get_Cext){
    mcSim->getCext();
    exit(0);
  }
  if(get_PeriodicPDB){
    mcSim->mainLattice->printPeriodicPDB();
    exit(0);
  }

  string fn_movie ="movie.pdb";
  if(recordMovie){
    moviestream.open(fn_movie.c_str());
  }
  Stats newStats;
  newStats.getLatticeStats(mcSim->mainLattice);
  cout<< "intial stats:"<<endl;
  newStats.printCout();
  //mainLattice->checkStats();

  //return(0);
}



