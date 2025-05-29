#include <math.h>

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
 
#include "EnergyMap.hpp"
#include "Stats.hpp"
#include "Lattice.hpp"
#include "Native.hpp"
#include "ptemp.h"
#include "App.hpp"
//#define AGGR   


EnergyMap::EnergyMap(Lattice * l1){
  l=l1;
  for(int i=0;i<MAX_PROCS;i++){
    for(int j=0;j<QMAX;j++){ 
      energyMap[i][j]=0;
    }
  }
}

void EnergyMap::setBetaID(int bID){
  betaID=bID;
}

int lowestEfound=8000;
void EnergyMap::mapStats(Stats & stats, double boltz){
  if(!App::nativeStats){
    if(stats.Cext < 0 || stats.Cext >=QMAX){
      cout<<"ERROR too many (" << stats.Cext << ") external contacts"<<endl;
      exit(1);
    }
    if(isnan(boltz)){
      cout<<"ERROR boltz nan"<<endl;
      exit(1);
    }
    energyMap[betaID][stats.Cext]+= boltz;
  }else{
    //check if within limits
    if(stats.Nint < 0 || stats.Nint >=QMAX){
      cout<<"ERROR too many (" << stats.Nint << ") native contacts"<<endl;
      exit(1);
    }
    // check if we need to find a lowest energy structure
    if(App::checkBelowNative && l->native->isBelowNativeEnergy(stats.Etot)){
      checkBelowNative(stats, boltz);
    }
    // map stats, depending on which type is requested
    if(App::cintStats){
      energyMap[betaID][stats.Cint] += boltz;
    }else{
      energyMap[betaID][stats.Nint] += boltz;
    }
  }
  // if collecting hbond stats, add these
  if(App::hbondStats){
    hbondMap[betaID][stats.Hext] +=boltz;
  }
  //add energy (for calculation CV)
  eMap[betaID][stats.Etot] +=boltz;
}

void EnergyMap::initTableStats(){
  std::stringstream fnStats;
  fnStats  << "statsTable"<<myRank <<".txt"; 
  tableStatsStream.open(fnStats.str().c_str());
  
  tableStatsStream<<"#step"<<"\t";
  tableStatsStream<<"T"<<"\t";

  tableStatsStream<<"beta"<<"\t";
  // tableStatsStream <<"boltz"<<"\t";
  
  tableStatsStream<<"weight"<<"\t";

  tableStatsStream <<"Eint"<<"\t";
  tableStatsStream <<"Eext"<<"\t";
  tableStatsStream <<"Esol"<<"\t";
  tableStatsStream <<"Etot"<<"\t";

  tableStatsStream <<"Cint"<<"\t";
  tableStatsStream <<"Cext"<<"\t";
  tableStatsStream <<"Ctot"<<"\t";

 if(App::nativeStats){
    tableStatsStream <<"Nint"<<"\t";
    tableStatsStream <<"Next"<<"\t";
    tableStatsStream <<"Ntot"<<"\t";
  }
  
  if(App::hbondStats){
    tableStatsStream <<"Hint"<<"\t";
    tableStatsStream <<"Hext"<<"\t";
    tableStatsStream <<"Htot"<<"\t";
  }

  if(App::eBetaMu>0 || App::clustStats){
    tableStatsStream <<"eBM"<<"\t";
    tableStatsStream <<"clust"<<"\t";
    tableStatsStream <<"maxCl"<<"\t";
    tableStatsStream <<"nChains"<<"\t";
  }

  tableStatsStream <<endl;
}

void EnergyMap::writeTableStats(Stats & stats,int step,double weight, ClusterInfo ci){
  double beta = betas[betaID]*100;
  std::stringstream buf;
  buf << std::fixed << std::setprecision(4) << 1/beta;  

  tableStatsStream<<step<<"\t";
  tableStatsStream<<buf.str()<<"\t";
  
  tableStatsStream <<beta <<"\t";
  //tableStatsStream <<boltz<<"\t";

  tableStatsStream <<weight <<"\t";

  tableStatsStream <<stats.Eint<<"\t";
  tableStatsStream <<stats.Eext<<"\t";
  tableStatsStream <<stats.Esol<<"\t";
  tableStatsStream <<stats.Etot<<"\t";

  tableStatsStream <<stats.Cint<<"\t";
  tableStatsStream <<stats.Cext<<"\t";
  tableStatsStream <<stats.Ctot<<"\t";

  if(App::nativeStats){
    tableStatsStream <<stats.Nint<<"\t";
    tableStatsStream <<stats.Next<<"\t";
    tableStatsStream <<stats.Ntot<<"\t";
  }
  
  if(App::hbondStats){
    tableStatsStream <<stats.Hint<<"\t";
    tableStatsStream <<stats.Hext<<"\t";
    tableStatsStream <<stats.Htot<<"\t";
  }
  if(App::eBetaMu>0 || App::clustStats){
    tableStatsStream <<App::eBetaMu<<"\t";
    tableStatsStream <<ci.clusters<<"\t";
    tableStatsStream <<ci.largestCluster<<"\t";
    tableStatsStream <<ci.totalChains<<"\t"; 
  }

  tableStatsStream <<endl;
}
void EnergyMap::finalizeTableStats(){
  tableStatsStream.close();
}


void EnergyMap::checkBelowNative(Stats & stats, double boltz){
  if(App::findLowestE){
    if(lowestEfound > stats.Etot){
      cout<< "energy below native energy"<<endl;
      cout<< "Mapped stats:"<<endl;
      stats.printCout();
      cout<< "Lattice stats:"<<endl;
      l->stats.printCout();
      if(l->stats.Etot == stats.Etot){
	std::stringstream outPDB;
	outPDB << "lowestEnergyN"<<myRank <<".pdb"; 
	l->writePDBMultiChain(outPDB.str());
	lowestEfound=stats.Etot;
	std::stringstream outSTATS;
	outSTATS  << "lowestEnergyN"<<myRank <<".txt"; 
	ofstream infoFile;
	infoFile.open(outSTATS.str().c_str());  
	infoFile<<stats.print2string()<<endl; 
	infoFile << "beta " << l->getBetaMoves() <<endl;
	infoFile.close();
      }
    }
  }else{
    if(! l->native->hasNativeStructure(stats.Ntot) ){ 
      cout<<"ERROR energy below native energy"<<endl;
      cout<< "Mapped stats:"<<endl;
      stats.printCout();
      cout<< "Lattice stats:"<<endl;
      l->stats.printCout();
      if(l->stats.Etot == stats.Etot){
	l->writePDBMultiChain("belowNativeEnergy.pdb");
	exit(1);
      }   
    }
  }
}

void EnergyMap::printGnuPlotDataWR(string fn_stats,int num_procs, int myRank){
  ofstream energyFile;
  energyFile.open(fn_stats.c_str());
  cout << "printing to "<<fn_stats <<endl;

  for(int betaID=0;betaID<num_procs;betaID++){ 
    energyFile<<"BETA\t"<<100.0 *  betas[betaID] << endl;
    
    for(int q=0;q<QMAX;q++){
      if(energyMap[betaID][q]>0){
	energyFile<< q <<"\t"<< energyMap[betaID][q] <<endl;
      }
    }
    //energyFile.close();
  }
  energyFile.close();
  cout << "finished printing to "<<fn_stats <<endl;
}




void EnergyMap::printEMap(string fn_stats,int num_procs, int myRank){
  
  ofstream energyFile;
  energyFile.open(fn_stats.c_str());
  cout << "printing to "<<fn_stats <<endl;
  for(int betaID=0;betaID<num_procs;betaID++){
    energyFile<<"BETA\t"<<100.0 *  betas[betaID] << endl;
    //__gnu_cxx::hash_map<int,double>::iterator iter;
    tr1::unordered_map<int,double>::iterator iter;
    for( iter = eMap[betaID].begin(); iter != eMap[betaID].end(); iter++ ) {
      energyFile<< (double) iter->first / 100.0 <<"\t"<< iter->second<<endl;
    }
  }
  energyFile.close();
}



void EnergyMap::printHBondMap(string fn_stats,int num_procs, int myRank){
  ofstream energyFile;
  energyFile.open(fn_stats.c_str());
  cout << "printing to "<<fn_stats <<endl;

  for(int betaID=0;betaID<num_procs;betaID++){ 
    energyFile<<"BETA\t"<<100.0 *  betas[betaID] << endl;
    
    for(int q=0;q<QMAX;q++){
      if(hbondMap[betaID][q]>0){
	energyFile<< q <<"\t"<< hbondMap[betaID][q] <<endl;
      }
    }
    //energyFile.close();
  }
  energyFile.close();
  cout << "finished printing to "<<fn_stats <<endl;
}

