#include "GranCan.hpp"
#include "MonteCarlo.hpp"
#include "Lattice.hpp"
#include "Chain.hpp"
#include "MolBox.hpp"
#include "EnergyMap.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "Cluster.hpp"
#include "ptemp.h"
#include "Moves.hpp"

#include "App.hpp"


#include <string.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

#define CORRECTION_0CHAINS


extern Pos local[6];

 

MonteCarlo::MonteCarlo(){
  mainLattice = new Lattice;
  //molbox=new MolBox(mainLattice);  // needed to for Stats::setLattice
  grancan = new GranCan(mainLattice);


  pGlobal =0.1;
  pInsDel =0.01;
  pClusterMove=0.005;

  pChangeStrandCoil=0.5;
  
  steps =10000;
  total_swaps=1;


}



int MonteCarlo::getCext(){
  mainLattice->stats.printCout();
  cout<<"CEXT "<<App::fn_in<<" "<<mainLattice->stats.getCext()<<endl;
  return mainLattice->stats.getCext();

}


int MonteCarlo::MC(double beta, int betaId, bool betaChanged){
  
  grancan->totalEmptySteps=0;
  //  cout<< "starting MC "<<endl;
  
  //mainLattice->writePDBMultiChain("lastConfig.pdb");
  mainLattice->setBetaMoves(beta);
  // cout<< "betas mainLattice set "<<endl;
  if(App::eBetaMu>0){
    grancan->molbox->createNewSet(beta,App::num_configs,betaChanged);
    //cout<< "new config  set created"<<endl;
  }
  mainLattice->energyMap->setBetaID(betaId);

  for(long i=0;i<steps;i++){
    //    cout<< i <<endl;
    //mainLattice->stats.printCout();
    if(drand48() < pGlobal){
      Chain * c= mainLattice->chains[(int) floor(mainLattice->nChains * drand48())];
      Moves::globalMove(mainLattice,c);
      //cout<<"finished global move"<<endl;
       // mainLattice->checkStats();
    }
    if(drand48() < pInsDel){   
      if(drand48()<0.5){
		grancan->trialInsertChain();
      }else{
		grancan->trialDeleteChain();
      }
    }
    
    if(drand48() < pChangeStrandCoil){
      Chain * c= mainLattice->chains[(int) floor(mainLattice->nChains * drand48())];
      Moves::changeStrandCoil(mainLattice,c);
      //cout<<"finished strand coil switch"<<endl;
      //mainLattice->checkStats();
    }else{
      Chain * c= mainLattice->chains[(int) floor(mainLattice->nChains * drand48())];
      Moves::shuffleSpin(mainLattice,c);
      //cout<<"finished spin flip"<<endl;
      //mainLattice->checkStats();
    }
    Chain * c= mainLattice->chains[(int) floor(mainLattice->nChains * drand48())];
    Moves::localMove(mainLattice,c);
     //cout<<"finished local move"<<endl;
     //mainLattice->checkStats();
  }

  if(App::clusterMoves && pClusterMove < drand48()){
    grancan->clusterMove();
  }


  for(int c=0;c<mainLattice->nChains;c++){
     mainLattice->chains[c]->checkAllFwdBkwd();
  }
  ClusterInfo ci;
  if(pInsDel >0 || App::clustStats){
    Cluster cl;
    ci = cl.checkCluster(mainLattice,beta);
  }
  mainLattice->checkStats();
  if(App::recordMovie){
    mainLattice->writePDBMultiChain(App::moviestream,App::moviestep);
    App::moviestep++;
  } 


  mainLattice->energyMap->writeTableStats(mainLattice->stats,App::tableStep,1.0,ci);
  if(grancan->totalEmptySteps>0){
    Stats emptyStats;
    double weight = grancan->totalEmptySteps/(double)steps;
    ClusterInfo ci_empty = {0,0,0};
    mainLattice->energyMap->writeTableStats(emptyStats,- App::tableStep,weight,ci_empty);
  }
  App::tableStep++;
  //mainLattice->stats.printCout();
  return mainLattice->stats.getEtot(); 
}



void MonteCarlo::getMCStats(int num_procs, int myRank){
  Cluster cl;
  cl.startCounting(mainLattice);
  cl.printCout();

  std::stringstream outPDB;
  outPDB << "outN"<<myRank<<".pdb";
  mainLattice->writePDBMultiChain(outPDB.str());


  std::stringstream outSTATS;
  outSTATS << App::fn_stats << "N"<<myRank <<".tab"; 
  mainLattice->energyMap->printGnuPlotDataWR(outSTATS.str(),num_procs,myRank);
  
  
  std::stringstream outEMAP;
  outEMAP << App::fn_eMap << "N"<<myRank <<".tab"; 
  mainLattice->energyMap->printEMap(outEMAP.str(),num_procs,myRank);

  if(App::hbondStats){
    std::stringstream outHBondMap;
    outHBondMap << App::fn_hbond << "N"<<myRank <<".tab"; 
    mainLattice->energyMap->printHBondMap(outHBondMap.str(),num_procs,myRank);
  }

  mainLattice->stats.printCout();
  cout<<"free chains: "<<grancan->getNumFreeChains()<<endl;

  mainLattice->energyMap->finalizeTableStats();
    
  //printGnuPlotData(fn_stats,num_procs,myRank);
  //printEMap(fn_eMap,num_procs,myRank);
  //transitionRates->print2file("transitions",myRank);
  //sqrg->print2file("sqrg",myRank);
}




