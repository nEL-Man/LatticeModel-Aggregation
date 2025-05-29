/**
 * \file   Design.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Design
 */


#include "Design.hpp"
#include "Lattice.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "AA.hpp"
#include "StatsMethods.hpp"
#include "App.hpp"
#include "AADistr.hpp"
#include "Moves.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>


double countAA[MAXAA];

double pChangeStrandCoil2=0.5;

double systemVar=0;

int numRes=0;


Design::Design(){
  mainLattice = new Lattice;
  mainLattice->setAA(App::fn_aa,App::indx_water); 
  //// set AA distribution
  if(App::AAdistribution){
    aaDistr = new AADistr();
    aaDistr->init(App::fn_AAdistribution,mainLattice->aaInt);
  }
}


void Design::init(string fn_in){
  mainLattice->readPDBMultiChain(fn_in);
  if(App::AAdistribution){
    aaDistr->setSequence(mainLattice);
  }else{
    initCountAA();
  }
  // mainLattice->stats.getLatticeStats(mainLattice);
  mainLattice->stats.printCout();

  /// set random spins
  for(int cn =0;cn<mainLattice->nChains;cn++){
    if(! mainLattice->chains[cn]->frozen && App::flipSpins)
      mainLattice->chains[cn]->setRandomSpins();
  }
}

void Design::initOldFormat(string fn_in){
  mainLattice->oldReadPDBMultiChain(fn_in);
  initCountAA();
  


  /// set random spins
  for(int cn =0;cn<mainLattice->nChains;cn++){
    if(! mainLattice->chains[cn]->frozen && App::flipSpins)
      mainLattice->chains[cn]->setRandomSpins();
  }
}



void Design::designProcedure(double beta, long numSteps){

  mainLattice->setBetaMoves(beta);
  cout<< "energy Beta "<<beta<<endl;

  //for(int cn =0;cn<mainLattice->nChains;cn++){
  //  if(! mainLattice->chains[cn]->frozen && App::flipSpins)
  //  mainLattice->chains[cn]->setRandomSpins();
  // }
  mainLattice->stats.getLatticeStats(mainLattice);

  for(long i=0;i<numSteps;i++){
    // change AA
    int cn = (int) floor(mainLattice->nChains * drand48());
    if(!mainLattice->chains[cn]->frozen){ 
      changeAA(cn);
      //if(i < (numSteps/2)){
      if(App::flipSpins){
	if (drand48() < pChangeStrandCoil2){
	  // change state
	  Chain * c= mainLattice->chains[(int) floor(mainLattice->nChains * drand48())];
	  Moves::changeStrandCoil(mainLattice,c);
	}else{ 
	  // change spin
	  Chain * c= mainLattice->chains[(int) floor(mainLattice->nChains * drand48())];
	  Moves::shuffleSpin(mainLattice,c);
	}
      }
    }
  }
  mainLattice->stats.printCout();
  mainLattice->checkStats();
  if(App::AAdistribution){
    aaDistr->check();
  }else{
    checkVariance();
  }
}


void Design::finalize(){
  std::stringstream fn_out_design_out;
  std::stringstream fn_out_design_pdb;
  
  fn_out_design_out << "design" <<myRank <<".out"; 
  fn_out_design_pdb << "design" <<myRank <<".pdb"; 
  
  writeDesignStats(fn_out_design_out.str());
  mainLattice->writePDBMultiChain(fn_out_design_pdb.str());
  
  mainLattice->stats.printCout();
  mainLattice->checkStats();
  if(App::AAdistribution){
    aaDistr->check();
  }else{
    checkVariance();
  }

}


void Design::checkVariance(){
 if(systemVar != getSystemVar()){
    cout <<"ERROR in design process accummulative variance:"<<systemVar;
    cout <<" does not match calc variance:" <<getSystemVar()<<endl;
    //exit(1);
 }
}

void Design::changeAA(int chain){
  Chain * c = mainLattice->chains[chain];
  int n = (int) floor(c->N * drand48());
  Residue * res = c->residues[n];
  int oldAA = res->aa;
  int newAA;
  do{  newAA = (int) floor(AA::NUMAA * drand48());}while(newAA == AA::designWater);
  
 
  //// draw AA acid with respect to correct distribution
  double localVar =0; double dist=0;
  if(App::AAdistribution){
    dist = aaDistr->getDistanceSubst(oldAA,newAA);
    double acc = exp(-dist / App::designT);
    if(!(drand48()<acc)){
      return;
    }
    
  }else{
    localVar = countAA[oldAA]/(countAA[newAA]+1.0);
    double accVar =pow(localVar ,App::designT);
    if(!(drand48()<accVar)){
      return;
    }
  }
  
  Stats oldLocalStats; 
  Stats newLocalStats; 
  oldLocalStats.localStats(res,res->pos,res->spin,res->state,oldAA,mainLattice);
  newLocalStats.localStats(res,res->pos,res->spin,res->state,newAA,mainLattice);
  Stats newLatticeStats =  mainLattice->stats.delta(newLocalStats,oldLocalStats);
  
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )*  mainLattice->getBetaMoves());
  bool accept=  dE<=0 || (drand48() < boltz);    
  
  if(accept){
    // update AA
    res->aa=newAA;
    // update energy
    mainLattice->stats = newLatticeStats;
    // update AAcount

    if(App::AAdistribution){
      aaDistr->substitute(oldAA,newAA,dist);
    }else{
      countAA[oldAA]--;
      countAA[newAA]++;
      systemVar = systemVar *localVar;
    }
  }
}



void Design::initCountAA(){
  numRes=0;
  for(int i=0;i<AA::NUMAA;i++){
    countAA[i]=0;
  }

  for(int nc=0;nc<mainLattice->nChains;nc++){
    Chain * ch = mainLattice->chains[nc];
    for (int n=0;n<ch->N;n++){
      Residue * res =ch->residues[n];
      countAA[res->aa]++;
      numRes++;
    }
  }
  systemVar = getSystemVar();
}


double factorial(double number) {
  double temp;
  if(number <= 1.0) return 1.0;
  temp = number * factorial(number - 1.0);
  return temp;
}

double Design::getSystemVar(){
  double ni=1.0;

  for(int aa=0;aa<AA::NUMAA;aa++){
    if(countAA[aa]){
      ni *= factorial(countAA[aa]);
    }
  }
  //cout<<usedRes<<endl;

  //return (factorial((double)numRes))/((double) ni);
  return (1.0)/((double) ni);
}


void Design::writeDesignStats(string s){
  ofstream outs;
  cout<< "printing to "<< s<<endl;
  outs.open(s.c_str()); 

  //designT
  outs<<"designT:"<<App::designT<<endl;
  //energyT
  outs<<"energyT:"<<mainLattice->getBetaMoves()<<endl;
  //sequence
  for(int nc=0;nc<mainLattice->nChains;nc++){
    Chain * ch = mainLattice->chains[nc];
    outs<<"chain no:"<< nc<<endl;
    for (int n=0;n<ch->N;n++){
      Residue * res =ch->residues[n];
      outs<<AA::int2aa[res->aa]<<" ";
    }
    outs<<endl;
  }
  //Stats internal energy 
  outs<< mainLattice->stats.print2string();
 
  //sequence variance
  if(App::AAdistribution){
    outs<< aaDistr->toString() <<endl;
  }else{
    outs<<"systemsVar:"<<systemVar<<endl;
    //AA count
    outs<<"AA count:"<< endl;
    for(int aa=0;aa<AA::NUMAA;aa++){
      outs <<AA::int2aa[aa]<<" "<<countAA[aa]<<", ";
    }
    outs<<endl;
  }
 outs.close();
}
