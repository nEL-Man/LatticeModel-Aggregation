
#include "GranCan.hpp"
#include "Lattice.hpp"
#include "Chain.hpp"
#include "MolBox.hpp"
#include "EnergyMap.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "Cluster.hpp"
#include "ptemp.h"

#include "App.hpp"





#include <string.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>




extern Pos local[6];

 

GranCan::GranCan(Lattice * ml){

  mainLattice = ml;
  molbox=new MolBox(mainLattice);  // needed to for Stats::setLattice

  Volume = (double) (LX * LY * LZ) ; 
  App::eBetaMu = -0.0000005;
}



int GranCan::trialInsertChain(){
  //cout<< "trial insert chain"<<endl;
  // find free chains
  int l=0;
  for(int c=0;c<mainLattice->nChains;c++){
    if(! mainLattice->chains[c]->hasCext() ){
      l++;
    }
  }  
  int freeChains=l;
 // accept depending on mu


  double prefactor;
  if(LAMBDA_COR){
    prefactor = 1.0/pow(mainLattice->getRealBeta(),3.0/2.0);
  }else{
    prefactor =1.0;
  }
  
  double accRatio =  prefactor*App::eBetaMu * (double) Volume/ ((double)freeChains +1.0);
  if(drand48() < accRatio){

    // get configuration
    Config * config =  molbox->getNewConfig();
    Config * newConfig= new Config;
    newConfig->nTotal = config->nTotal;
    // add random translation to points (could do this in molbox) 
    Pos addP;
    addP.x = (int) floor( LX * drand48());
    addP.y = (int) floor( LY * drand48());
    addP.z = (int) floor( LZ * drand48()); 
    
    for(int n=0;n<config->nTotal;n++){
      //cout<<"config res:"<<n<<endl;
      Pos posN = config->positions[n]+addP;
      posN.periodicBoundary();
      // check if lattice points not occupied and not bound
      if(mainLattice->getResidue(posN) != NULL){
	delete newConfig;
	return -1;
      }
      for(int k=0;k<6;k++){
	Pos posNB = local[k] + posN;
	posNB.periodicBoundary();
	// check if bound, threonine (aa==13) contact does not count
	//if(mainLattice->getResidue(posNB)!= NULL && mainLattice->getResidue(posNB)->aa != 13 ) {
	if(mainLattice->getResidue(posNB)){
	  delete newConfig;
	  return -1;
	}
      }
#ifdef CROSSNBS
      Pos spin = config->spins[n];
      Pos * crossNbs = Stats::getCrossNbsSpinDir(spin);
      for(int k=0;k<4;k++){
	Pos posNB =crossNbs[k]+ posN;
	posNB.periodicBoundary();
	Residue * resNB = mainLattice->getResidue(posNB);
	if(resNB!=NULL){
	  bool sameChain = false;
	  if(!sameChain){
	    bool spinsOpposite= (resNB->spin == (Pos(0,0,0)- spin));
	    if(spinsOpposite){
	      delete newConfig;
	      return -1;
	    }
	  }
	}
      }
      //delete crossNbs;
#endif
      //Note this could be more efficient
      newConfig->positions[n]= posN;
      newConfig->spins[n]=config->spins[n];
      newConfig->states[n]=config->states[n];
    }
    
    
    // insert, by creating new chain
    mainLattice->insertChain(newConfig,molbox->sequence,mainLattice->getBetaMoves());
    delete newConfig;
    // mainLattice->writePDBMultiChain("newlyinsertedchain.pdb");
    return 1;
  }
  return 0;
}

// this may become slow for fully aggregated box
int GranCan::trialDeleteChain(){
  // cout<< "trial delete chain"<<endl;
  int freeChns[MAX_CHAINS];
  // find free chains
  int l=0;
  for(int c=0;c< mainLattice->nChains;c++){
    if(! mainLattice->chains[c]->hasCext()){
      freeChns[l]= c;
      l++;
    }
  }  
  int freeChains=l;
  if(freeChains < 1) return -2;
  if(mainLattice->nChains <=1 ){ 
#ifdef CORRECTION_0CHAINS
    double prefactor;
    if(LAMBDA_COR){
       prefactor = pow(mainLattice->getRealBeta(),3.0/2.0);
    }else{
      prefactor =1.0;
    }

   
    double accRatio = prefactor*((double)freeChains ) / (App::eBetaMu *  Volume);
    if(drand48() < accRatio){
      
      double alpha = App::eBetaMu * ((double) Volume)/prefactor;
      //double r = drand48();
      double num_steps;
      if(alpha > 0.9999){
	num_steps=1;
      }else{
	num_steps = ceil(log(drand48())/log(1.0-alpha));
      }
      double moves = (1.0+ pGlobal)*(2*num_steps/pInsDel);
      totalEmptySteps += moves;
      //cout<<"Volume "<<Volume<<endl;
      //cout<<"r "<<r<<endl;
      //cout<<"alpha "<<alpha<<endl;
      //cout<<"log(1.0-alpha)"<<log(1.0-alpha)<<endl;
      //cout<<"pG "<<pGlobal<<endl;
      //cout<<"n_s "<<num_steps<<endl;
      //cout<<"moves "<<moves<<endl;
      Stats emptyStats ;
      mainLattice->energyMap->mapStats(emptyStats,moves);
      
    }
    return -1;
#else
    return -1;
#endif
  }
  /*do{  
    chain = (int) floor(drand48() * nChains);  
    }while(chains[chain]->Cext!=0);
  */
	
  // accept deletion depending on mu
  double prefactor;
  if(LAMBDA_COR){
    prefactor = pow(mainLattice->getRealBeta(),3.0/2.0);
  }else{
    prefactor =1.0;
  }


  double accRatio = prefactor*((double)freeChains ) / (App::eBetaMu *  Volume);
  if(drand48() < accRatio){
     // pick a free chain at random
   
    int chain= freeChns[(int) floor(drand48() * freeChains)];
    // delete chain from lattice
    mainLattice->deleteChain(chain);
    return 1;
  }
  return 0;
}

int GranCan::getNumFreeChains(){
  // cout<<"getNumfree chains" <<endl;
 int l=0;
  for(int c=0;c< mainLattice->nChains;c++){
    if(! mainLattice->chains[c]->hasCext()){
      l++;
    }
  } 
  return l;
}





int GranCan::clusterMove(){
  /// calculate clusters
  Cluster cl;
  cl.startCounting(mainLattice);
  int numClusts = cl.getTotalClusters();
  //  cl.printCout();

  /// pick random clusters
  int cluster2move = (int) floor(numClusts * drand48());

  /// get all chains in cluster
  vector<int> chainList = cl.getChains(cluster2move);
  int numChains = chainList.size();

  /// set rotation and translation
  //int dir_translation= (int)floor(6 * drand48());
  //int magnitude_translation = (int)floor(MAX_TRANSLATION * drand48());
  int rotationalDir= (int)floor(6 * drand48());

  /// choose random pivot on lattice
  int xpivot =  (int)floor(LX * drand48()); 
  int ypivot =  (int)floor(LY * drand48());
  int zpivot =  (int)floor(LZ * drand48());
  Pos posPivot =  Pos(xpivot,ypivot,zpivot);

  vector< vector <Pos> > newPos (numChains, vector<Pos>());

  /// for each chain, perform rotation and translation,
  /// while keeping old coordnates
  /// check for clashes or touches: reject, restore old coordinates and sample
  for(unsigned int i=0; i< chainList.size();i++){
    int cn = chainList[i];
    Chain * ch = mainLattice->chains[cn];
    for(int n=0; n<ch->N; n++){
      Residue * res = ch->residues[n];
      Pos pp = ch->rotationPosition(res->pos,rotationalDir,posPivot);
      newPos[i].push_back(pp);
      if(mainLattice->checkForClashesAndTouches(newPos[i][n])){
	/// reject move
	// l->sampleCurrentStats();
	return -1;
      }      

    }
  }
  /// now we know there are no clashes our touches
  /// and the move can be performed
  /// note that the move will alsways get accepted, since dE=0

  /// Perform move, update lattice
  for(unsigned int i=0; i< chainList.size();i++){
    int cn = chainList[i];
    Chain * ch = mainLattice->chains[cn];
    for(int n=0; n<ch->N;n++){
      Residue * res = ch->residues[n];
      Pos newSpin = ch->rotationSpin(res->spin,rotationalDir);
      mainLattice->emptyLatticePos(res->pos);

      res->pos = newPos[i][n];
      res->spin = newSpin;
      mainLattice->setResidue(newPos[i][n],res);
    }
    ch->setFwdBkwdSpins();
  }
  return 1;
}


