#include "MolBox.hpp"

#include "Chain.hpp"
#include "Lattice.hpp"
#include "Native.hpp"
#include "Stats.hpp"
#include "App.hpp"
#include "Moves.hpp"

MolBox::MolBox(Lattice * mainL){
  //beta = 0.08;
  currentConfig =0;
  numConfigs =1000;
  box = new Lattice;
  mainLattice = mainL;
  for(int i =0;i<MAX_CONFIGS;i++){
    configSet[i]=NULL;
  }
}

void MolBox::setAA(string fn, int indx){
  box->setAA(fn,indx);
}

/*void MolBox::setBeta(double b){
  beta=b;
  }*/

double MolBox::getBeta(){
  return box->getBetaMoves();
}


void MolBox::setMolecule(Chain * chain){
  // lattice should be clean 

  molecule = new Chain(box,0);
  box->chains[0]=molecule;
  molecule->frozen = false;
  molecule->locked = false;
  molecule->N = chain->N;
  
  // copy and insert each residue
  for(int n=0;n<chain->N;n++){
    Residue * res = new Residue;
    res->pos.x = chain->residues[n]->pos.x;
    res->pos.y =chain->residues[n]->pos.y;
    res->pos.z =chain->residues[n]->pos.z;
    res->spin = chain->residues[n]->spin;
    res->state = chain->residues[n]->state;
    res->aa =chain->residues[n]->aa;
    sequence[n]=res->aa;
    res->n =chain->residues[n]->n;
    res->chainNum=0;  
    molecule->residues[n]=res;
    box->r[res->pos.x][res->pos.y][res->pos.z] = molecule->residues[n];
  }
  molecule->setFwdBkwdSpins();
  box->nChains = 1;
  if(App::nativeStats){
    box->native = new Native(box);
  }
  box->stats.getLatticeStats(box);
  box->energyMap = NULL;// should not record eMap
}



void MolBox::createNewSet(double b,int num_configs, bool betaChanged){
  // need to set Stats::l
  //cout<< "CREATING NEW SET "<<endl;
  box->setBetaMoves(b);
  //molecule->setBeta(beta);

  if(betaChanged){   
    //cout<<"EQUIL BETA CHANGE"<<endl;
    equilibrate(EQUIL_STEPS);
  }
 
  currentConfig=0;
  if(num_configs > MAX_CONFIGS)num_configs = MAX_CONFIGS;
  numConfigs=num_configs;
 
  for(int nc =0;nc<numConfigs;nc++){
    //cout<< nc <<endl;
    equilibrate((int) floor(100 * drand48()));
    recordConfig(nc);         
  } 
}

void  MolBox::printConfig(int c){
  Config * config = configSet[c];
  cout<<"configuration"<< c<<endl;
  for(int n=0;n<config->nTotal;n++){
    cout<<"config res:"<<n<<" "<<config->positions[n].x<<" "<<config->positions[n].y<<" "<<config->positions[n].z<<endl;
  }
}



Config * MolBox::getNewConfig(){
  currentConfig++;
  if(currentConfig >= numConfigs){
    createNewSet(box->getBetaMoves(),numConfigs, false);
  }
  // perhaps better to choose @random?
  int c = (int) floor((numConfigs-1) * drand48());
  //cout<<"config no: "<<c<<endl;
  //printConfig(c);
  return configSet[c ];  
}
 

void MolBox::equilibrate(int nsteps){
  //cout<<"EQUILIBIRATING"<<endl;
  for(int i =0;i<nsteps;i++){
    //molecule->localMove();
    Moves::localMove(box,molecule);
    if(drand48() < 0.05 ){
      Moves::rotatePoint(box,molecule);
    }
    if(drand48() < 0.5){
      Moves::changeStrandCoil(box,molecule);
    }else{
      Moves::shuffleSpin(box,molecule);
    }
  }
}


void MolBox::recordConfig(int insertAt){
  //cout<<"recording config"<<endl;
  if(configSet[insertAt] !=NULL){
    delete configSet[insertAt];
    configSet[insertAt] = NULL;
  }
  
  Config * cc = new Config;
 
  for(int n=0; n<molecule->N;n++){
    // perhaps needs normalising ?

    cc->positions[n] = molecule->residues[n]->pos;
    cc->spins[n] =  molecule->residues[n]->spin;
    cc->states[n] = molecule->residues[n]->state;
    cc->nTotal = molecule->N;
  }

  configSet[insertAt]= cc;
}

