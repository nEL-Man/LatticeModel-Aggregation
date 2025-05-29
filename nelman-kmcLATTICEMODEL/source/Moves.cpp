#include "Chain.hpp"
#include "Moves.hpp"
#include "Lattice.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "EnergyMap.hpp"

#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

/******************Standard Move **************************/
/* Generate new configuration                        */
/* Check if newConfig illegal,if illegal return -1   */
/* Calculate energy difference between states        */
/* Sample boths states for statistics                */
/* accept or reject move                             */
/* reject: return 0                                  */
/* accept: update lattice and chain,return 1         */
/**********************************************************/



// class variables to store changes in stats for speedup
// cleared for every move
Stats Moves::oldLocalStats;
Stats Moves::newLocalStats;
Stats Moves::newLatticeStats;




int Moves::shuffleSpin(Lattice *l,Chain * c){
  int n =(int) floor(c->N * drand48());
  Residue * res = c->residues[n];
  //cout<<"shuffle spin "<<n<<endl;
  if(res->state==inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  Pos newSpin = local[(int) floor(6 * drand48())];
  if(newSpin == res->bkwd || newSpin == res->fwd){
    l->sampleCurrentStats();
    return -1;
  }  				  
  //energy acceptance + sampling
  oldLocalStats.clean(); 
  newLocalStats.clean();
 
  oldLocalStats.localStats(res,res->pos,res->spin,res->state ,l);
  newLocalStats.localStats(res,res->pos,newSpin,res->state,l);
  
  newLatticeStats.clean();
  newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  double bb = boltz/(boltz+1.0);
  if(l->energyMap != NULL){
    l->energyMap->mapStats(l->stats,1.0 - bb);
    l->energyMap->mapStats(newLatticeStats,bb);
  }
 
  if(dE<=0 || drand48() < boltz){
    //accept
    //cout<<"change state accepted"<<endl;
    res->spin = newSpin;
    l->stats = newLatticeStats;
    return 1;
  }else{
    return 0;
  }
}


int Moves::changeStrandCoil(Lattice *l,Chain * c){
  int n = (int) floor(c->N * drand48());
  Residue *res=c->residues[n];
  //cout<<"changeStrandCoil spin "<<n<<endl;
  State newState;
  if(res->state ==inCoil){
    if(c->strandPossible(res)){
      newState=inStrand;
    }else{
      l->sampleCurrentStats();
      return -1;
    }
  }else{
    newState=inCoil;
  }

  oldLocalStats.clean(); 
  newLocalStats.clean(); 
  oldLocalStats.localStats(res,res->pos,res->spin,res->state ,l);
  newLocalStats.localStats(res,res->pos,res->spin,newState,l);
  
  newLatticeStats.clean();
  newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  double bb = boltz/(boltz+1.0);
  if(l->energyMap != NULL){
    l->energyMap->mapStats(l->stats,1.0 - bb);
    l->energyMap->mapStats(newLatticeStats,bb);
  }
 
  if(dE<=0 || drand48() < boltz){
    //accept
    res->state = newState;
    //cout<<"change state accepted"<<endl;
    l->stats = newLatticeStats;
    return 1;
  }else{
    return 0;
  }
}




int Moves::localMove(Lattice *l,Chain * c){
  Pos posNew;
  Pos newSpin;
  int px=-1,py=-1,pz=-1; 
  int n = (int) floor(c->N * drand48());
  if(c->residues[n]->state == inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  if(n>0 && c->residues[n-1]->state == inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  if(n<(c->N-1) && c->residues[n+1]->state == inStrand){
    l->sampleCurrentStats();    
    return -1;
  }
  Residue * res=c->residues[n];
 

  //cout <<"local "<<n<<endl;
  if(n<0 ||n >= c->N){ cout<<"residue not in range"<<n<<endl; exit(-1);} 
  if(n==0 || n == c->N - 1){ // end of chain choose one of 6 directions directions
    int nn =1;
    if(n==(c->N-1)){nn = c->N-2;} // senearest neighbour   
    int rr= (int)floor(6 * drand48()); 
    // OPT slow, Behnaz: put links instead
    posNew = local[rr] + c->residues[nn]->pos;
    posNew.periodicBoundary();
    if(res->pos != posNew){
      newSpin = c->newSpinEndMove(res->pos,posNew,c->residues[nn]->pos,res->spin);
    }else{
      newSpin=res->spin;
    }
  }else{ // try cornerflip
    px = (res->pos.x == c->residues[n-1]->pos.x) & (res->pos.x == c->residues[n+1]->pos.x);
    py = (res->pos.y == c->residues[n-1]->pos.y) & (res->pos.y == c->residues[n+1]->pos.y);
    pz = (res->pos.z == c->residues[n-1]->pos.z) & (res->pos.z == c->residues[n+1]->pos.z);
    if(px+py+pz==2){
      // on a straight line
      l->sampleCurrentStats();
      return -1;
    }else if(px+py+pz==1){ 
      if(px){
	posNew.x=res->pos.x;
	posNew.y =  c->residues[n-1]->pos.y + c->residues[n+1]->pos.y - res->pos.y;
	posNew.z =  c->residues[n-1]->pos.z + c->residues[n+1]->pos.z - res->pos.z;
      }else if(py){
	posNew.y=res->pos.y;
	posNew.x =  c->residues[n-1]->pos.x + c->residues[n+1]->pos.x - res->pos.x;
 	posNew.z =  c->residues[n-1]->pos.z + c->residues[n+1]->pos.z - res->pos.z;
      }else{
	posNew.z=res->pos.z;
	posNew.x =  c->residues[n-1]->pos.x + c->residues[n+1]->pos.x - res->pos.x;
	posNew.y =  c->residues[n-1]->pos.y + c->residues[n+1]->pos.y - res->pos.y;
      }
      if(res->pos != posNew){
	newSpin=c->newSpinPosCornerFlip(res->spin,c->residues[n-1]->pos,c->residues[n+1]->pos);
      }else{
	newSpin=res->spin;
      }
    }else{
      cout << "break in chain: "<< c->chainNum <<endl;;
      cout << "trying to move:"<<n<< " with coords:" ;
      cout <<res->pos.x<<" "<<res->pos.y<<" "<<res->pos.z<<" "<<endl;;
      cout <<c->residues[n-1]->pos.x<<" "<<c->residues[n-1]->pos.y<<" "<<c->residues[n-1]->pos.z<<" "<<endl;
      cout <<c->residues[n+1]->pos.x<<" "<<c->residues[n+1]->pos.y<<" "<<c->residues[n+1]->pos.z<<" "<<endl;
      //printChain();
      l->writePDBMultiChain("chainBreak.pdb");
      //return -1;
       exit(1);
    }


  }
  //check if no steric clash main chain
  if(l->getResidue(posNew)!=NULL){
    if(l->getResidue(posNew)->n == n-2 && l->getResidue(posNew)->chainNum == c->chainNum){
      return crankshaft(l,c,n-1,n,px,py,pz);
    }else if(l->getResidue(posNew)->n == n+2 && l->getResidue(posNew)->chainNum== c->chainNum){
      return crankshaft(l,c,n,n+1,px,py,pz);
    }else{
      l->sampleCurrentStats();
      return(-1);
    }
  }
  Pos newFwd,newFwdPrev;
  Pos newBkwd,newBkwdNxt;
  //check if no violation n-1 spin and n+1 spin
  if(n>0) {
    newBkwd.periodicSubtraction(c->residues[n-1]->pos,posNew);
    newFwdPrev= Pos(0,0,0)- newBkwd;
    if(c->residues[n-1]->spin == newFwdPrev){
      l->sampleCurrentStats();
      return -1;
    }
  }

  if(n< c->N - 1){
    newFwd.periodicSubtraction(c->residues[n+1]->pos,posNew);
    newBkwdNxt = Pos(0,0,0)-  newFwd;
    if(c->residues[n+1]->spin == newBkwdNxt){
      l->sampleCurrentStats();
      return -1;
    }
  }


  //update lattice
  l->emptyLatticePos(res->pos);
  
  oldLocalStats.clean(); 
  newLocalStats.clean(); 
  //calculate old stats
  oldLocalStats.localStats(res,res->pos,res->spin,res->state,l);
  oldLocalStats.solventStats(posNew,l);

  //calculate new stats
  newLocalStats.localStats(res,posNew,newSpin,newFwd,newBkwd,res->state,res->aa,l);
  newLocalStats.solventStats(res->pos,l);

  newLatticeStats.clean();
  newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  double bb = boltz/(boltz+1.0);
  if(l->energyMap != NULL){
    l->energyMap->mapStats(l->stats,1.0 - bb);
    l->energyMap->mapStats(newLatticeStats,bb);
  }
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  if(accept){
    
    // update spins
    if(n>0){
      res->bkwd=newBkwd;
      c->residues[n-1]->fwd=newFwdPrev;
    }
    if(n<c->N-1){
      res->fwd=newFwd;
      c->residues[n+1]->bkwd = newBkwdNxt;
    }
    res->spin=newSpin;
    
    // update position and lattice 
    res->pos = posNew;
    l->setResidue(posNew,res);
  
    // update lattice stats
    l->stats = newLatticeStats;
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    l->setResidue(res->pos,res);  
    return(0);
  }
}



int Moves::crankshaft(Lattice *l,Chain *c, int n1, int n2,int px, int py, int pz){
  //cout << "crankshaft "<<n1<<" "<<n2<<endl;
  int n0,n3;
  Residue * resn1 = c->residues[n1];
  Residue * resn2 = c->residues[n2];

  //copy old residues
  // Residue * newRes1 = new Residue(resn1);
  //Residue * newRes2 = new Residue(resn2);

  if(resn1->state == inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  if(resn2->state == inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  // TO DO rewrite this:
  if(n1<n2){
    n0=n1-1;
    n3=n1+2;
  }else{
    n0=n1+1;
    n3=n1-2;
    //swap all
    int tmp =n0;
    n0=n3;
    n3=tmp;
    tmp=n1;
    n1=n2;
    n2=tmp;
  }
  Residue * resn0 = c->residues[n0];
  Residue * resn3 = c->residues[n3];
  if(n0>0 && resn0->state == inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  if(n3 <(c->N-1) && resn3->state == inStrand){
    l->sampleCurrentStats();
    return -1;
  }
  if(n0<0 || n3<0 || n0>= c->N || n3>= c->N){
    l->sampleCurrentStats();
    return -3;
  }
  // OPT only need dirOrthogonal
  int dirStalk,dirLegs,dirOrthogonal;
  if(px){
    dirOrthogonal = XDIR;
    if((resn1->pos.y - resn2->pos.y)==0){
      dirStalk=YDIR;
      dirLegs=ZDIR;
    }else{
      dirStalk=ZDIR;
      dirLegs=YDIR;
    }
  }else if(py){
    dirOrthogonal = YDIR;
    if((resn1->pos.x - resn2->pos.x)==0){
      dirStalk=XDIR;
      dirLegs=ZDIR;
    }else{
      dirStalk=ZDIR;
      dirLegs=XDIR;
    }
  }else {
    dirOrthogonal = ZDIR;
    if((resn1->pos.x - resn2->pos.x)==0){
      dirStalk=XDIR;
      dirLegs=YDIR;
    }else{
      dirStalk=XDIR;
      dirLegs=YDIR;
    }
  }
  int dirFlip;
  if( drand48()<0.5){
    dirFlip=1;
  }else{
    dirFlip=-1;
  }
 
  Pos old0 =Pos(resn0->pos.x,resn0->pos.y,resn0->pos.z);
  Pos old3 =Pos(resn3->pos.x,resn3->pos.y,resn3->pos.z);

  Pos newPos1;
  Pos newPos2;

  newPos1[dirOrthogonal]= old0[dirOrthogonal] + dirFlip;
  newPos2[dirOrthogonal]= old3[dirOrthogonal] + dirFlip;

  newPos1[dirStalk]= old0[dirStalk];
  newPos2[dirStalk]= old3[dirStalk];
 
  newPos1[dirLegs]= old0[dirLegs];
  newPos2[dirLegs]= old3[dirLegs];

  newPos1.periodicBoundary();
  newPos2.periodicBoundary();

  

  // test for collision
  if(l->getResidue(newPos1)!=NULL || 
     l->getResidue(newPos2)!=NULL){
    l->sampleCurrentStats();
    return -1;
  }


  
  // test for collision new spins with new fwd and new bkw
  Pos newBkwdn1,newFwdn0;
  Pos newFwdn2,newBkwdn3;
  newBkwdn1.periodicSubtraction(resn0->pos, newPos1);
  newFwdn2.periodicSubtraction(resn3->pos,newPos2);
  newFwdn0 = Pos(0,0,0)- newBkwdn1;
  newBkwdn3= Pos(0,0,0)- newFwdn2;
  if(resn0->spin==newFwdn0){l->sampleCurrentStats(); return -1;}
  if(resn3->spin==newBkwdn3){l->sampleCurrentStats(); return -1;}
  Pos newSpin1 = c->newSpinEndMove(resn1->pos, newPos1, resn0->pos, resn1->spin);
  Pos newSpin2 = c->newSpinEndMove(resn2->pos, newPos2, resn3->pos, resn2->spin);;   
 
   // remove residues involved in move from lattice 
  l->emptyLatticePos(resn1->pos);
  l->emptyLatticePos(resn2->pos);
  
  //old Stats
  oldLocalStats.clean();
  
  oldLocalStats.localStats(resn1,l);
  oldLocalStats.localStats(resn2,l);

  oldLocalStats.solventStats(newPos1,l);
  oldLocalStats.solventStats(newPos2,l);


  //calculate new stats
  newLocalStats.clean();
  newLocalStats.localStats(resn1,newPos1,newSpin1,resn1->fwd,newBkwdn1,resn1->state,resn1->aa,l);
  newLocalStats.localStats(resn2,newPos2,newSpin2,newFwdn2,resn2->bkwd,resn2->state,resn2->aa,l);

  newLocalStats.solventStats(resn1->pos,l);
  newLocalStats.solventStats(resn2->pos,l);

  newLatticeStats.clean();
  newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  int dE= newLocalStats.getDeltaE(oldLocalStats);

  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  double bb = boltz/(boltz+1.0);
  if(l->energyMap != NULL){
    l->energyMap->mapStats(l->stats,1.0 - bb);
    l->energyMap->mapStats(newLatticeStats, bb);
  }
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  if(accept){
    //cout<<"crankshaft accepted"<<endl;
    //update spins
    resn1->spin=newSpin1;
    resn2->spin=newSpin2;

    resn0->fwd=newFwdn0;
    resn1->bkwd=newBkwdn1;
    resn2->fwd=newFwdn2;
    resn3->bkwd=newBkwdn3;

    // update position and lattice 
    resn1->pos = newPos1;
    resn2->pos = newPos2;
    l->setResidue(newPos1,resn1);
    l->setResidue(newPos2,resn2);

    // update lattice stats
    l->stats = newLatticeStats;
    
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    // reset lattice 
    l->setResidue(resn1->pos,resn1);
    l->setResidue(resn2->pos,resn2);

    return(0);
  }


}


int Moves::globalMove(Lattice *l,Chain *c){
  if(drand48() < 0.5){
    return translate(l,c);
  }else{
    return rotatePoint(l,c);
  }
}


int Moves::translate(Lattice *l,Chain * c){
  int dir= (int)floor(6 * drand48());

  Pos posNew[MAX_RES];
  
  for(int n=0; n< c->N;n++){
    Residue * res = c->residues[n];
    posNew[n] = res->pos+ local[dir];
    posNew[n].periodicBoundary();
    Residue * testRes= l->getResidue(posNew[n]);
    //check for collision
    if(testRes!=NULL && testRes->chainNum != res->chainNum ){
#ifdef DEBUG
      cout << "clash with ligand "<<r[xnew][ynew][znew]->n <<endl;;
#endif
      l->sampleCurrentStats();
      return -1;
    }  
  }

  
  oldLocalStats.clean();
  newLocalStats.clean();

  //clear all positions involved in move from lattice
  for(int n=0;n<c->N;n++){
    l->emptyLatticePos(c->residues[n]->pos);
  }
  
  
  //calculate old and new stats
  //calculate solvents stats of new for old and old for new, not including self
  for(int n=0; n<c->N;n++){
    oldLocalStats.localStats(c->residues[n],l);
    oldLocalStats.solventStats(posNew[n],l);
    newLocalStats.localStats(c->residues[n],posNew[n],l);
    newLocalStats.solventStats(c->residues[n]->pos,l);
  }
  
  newLatticeStats.clean();
  newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);
  int dE= newLocalStats.getDeltaE(oldLocalStats);
  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  double bb = boltz/(boltz+1.0);
  if(l->energyMap != NULL){
    l->energyMap->mapStats(l->stats,1.0 - bb);
    l->energyMap->mapStats(newLatticeStats, bb);
  }
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  
  if(accept){
    // set positions residues
    // and set new lattice positions
    for(int n=0;n<c->N;n++){
      c->residues[n]->pos=posNew[n];
      l->setResidue(posNew[n],c->residues[n]);
    }

    l->stats = newLatticeStats;
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    //reset lattice
    for(int n=0; n<c->N;n++){
      l->setResidue(c->residues[n]->pos,c->residues[n]);
    } 
    return(0);
  }

}


int Moves::rotatePoint(Lattice *l,Chain * c){
  //cout<<"rotate point"<<endl;
  int rotationalDir= (int)floor(6 * drand48());
  int pivot= (int)floor(c->N * drand48());
  if(c->residues[pivot]->state == inStrand){l->sampleCurrentStats(); return -1;}
  int startRes=0;
  int endRes=c->N-1;
  // OPT don't need to include pivot
  if(!c->frozen) {// if frozen only allow global moves
    if(drand48()<0.5){
      startRes= pivot;
      if(pivot < (c->N - 1) && c->residues[pivot+1]->state==inStrand){l->sampleCurrentStats(); return -1;}
    }else{
      endRes=pivot;
      if(pivot > 0 && c->residues[pivot-1]->state==inStrand){l->sampleCurrentStats(); return -1;}
    }
  }
  //cout<<"pivot"<<pivot<<endl;
  //cout<<"startRes"<<startRes<<endl;
  //cout<<"endRes"<<endRes<<endl;
  Pos posNew[MAX_RES];
  // test for collisions
  if(startRes==endRes){return -1;}

  for(int n=startRes; n<=endRes;n++){
    Residue * res = c->residues[n];
    posNew[n]= c->rotationPosition(res->pos,rotationalDir,c->residues[pivot]->pos);
#ifdef DEBUG
    cout<<"move to (" <<chainNew[n]->pos.x<<","<<chainNew[n]->pos.y<< ","<<chainNew[n]->pos.z<<")"<<endl;
#endif
    Residue * testRes= l->getResidue(posNew[n]);
    if(testRes!=NULL && (testRes->chainNum != res->chainNum || testRes->n<startRes || testRes->n> endRes  )){
      l->sampleCurrentStats();
      return -1;
    }  
  }
  //check for spin new fwd bkwd violation
  //TO DO perhap faster if above previous block
  Pos newFwdPivot;
  Pos newBkwdPivot;
  
  if(startRes==pivot){
    newFwdPivot.periodicSubtraction(posNew[pivot+1],posNew[pivot]);
    newBkwdPivot =c->residues[pivot]->bkwd;
    if(c->residues[pivot]->spin == newFwdPivot){l->sampleCurrentStats(); return -1;}
  }else{
    newBkwdPivot.periodicSubtraction(posNew[pivot-1],posNew[pivot]);
    newFwdPivot = c->residues[pivot]->fwd;
    if(c->residues[pivot]->spin == newBkwdPivot){l->sampleCurrentStats();return -1;}
  }
  

  
  //calculate new Spins
  Pos newSpins[MAX_RES];
  Pos newFwd[MAX_RES];
  Pos newBkwd[MAX_RES];
  for(int n=startRes; n<=endRes;n++){
    newSpins[n]= c->rotationSpin(c->residues[n]->spin,rotationalDir);
    newFwd[n] = c->rotationSpin(c->residues[n]->fwd,rotationalDir);
    newBkwd[n] =c->rotationSpin(c->residues[n]->bkwd,rotationalDir);
  }
  newSpins[pivot]=c->residues[pivot]->spin; // don't change spin pivot
  if(pivot != (c->N-1)){
    newFwd[pivot]= newFwdPivot;
  }else{
    newFwd[pivot]= Pos(0,0,0);
  }

  if(pivot !=0){
    newBkwd[pivot]= newBkwdPivot;
  }else{
    newBkwd[pivot]=Pos(0,0,0);
  }

  oldLocalStats.clean();
  newLocalStats.clean();
  
  int exclStart =startRes;
  int exclEnd = endRes;
  if(startRes==pivot){
    exclStart = startRes +1;
  }else{
    exclEnd = endRes -1;
  }

  // clear lattice positions involved in move
  for(int n=exclStart; n<=exclEnd;n++){
    l->emptyLatticePos(c->residues[n]->pos);
  }

  for(int n=exclStart; n<=exclEnd;n++){ 
    Residue * res = c->residues[n];
    oldLocalStats.localStats(res,res->pos,l);
    newLocalStats.localStats(res,posNew[n], newSpins[n],newFwd[n],newBkwd[n],res->state,res->aa,l);
    oldLocalStats.solventStats(posNew[n],l);
    newLocalStats.solventStats(res->pos,l);
  }

  newLatticeStats.clean();
  newLatticeStats =  l->stats.delta(newLocalStats,oldLocalStats);

  int dE= newLocalStats.getDeltaE(oldLocalStats);

  double boltz = exp(((-(float) (dE)) )* l->getBetaMoves());
  double bb = boltz/(boltz+1.0);
  if(l->energyMap != NULL){
    l->energyMap->mapStats(l->stats,1.0 - bb);
    l->energyMap->mapStats(newLatticeStats, bb);
  }
  bool accept=false;
  if(dE<=0 || drand48() < boltz){
    accept=true;
  }
  
  if(accept){
    //cout <<"ACCEPT"<<endl;
     
    for(int n=startRes;n<=endRes;n++){
      c->residues[n]->pos=posNew[n];
      l->setResidue(posNew[n],c->residues[n]);
    }  
    
    // set spins and fwd bkwd  
    for(int n=startRes;n<=endRes;n++){
      c->residues[n]->spin = newSpins[n];
      c->residues[n]->fwd = newFwd[n];
      c->residues[n]->bkwd = newBkwd[n];
      //setFwdBkwd(c->residues[n]);
    }
    
    l->stats = newLatticeStats;
    
    if(dE==0){
      return(1);
    }else{
      return(2);
    }
  }else{
    for(int n=startRes;n<=endRes;n++){
      //reset lattice
      l->setResidue(c->residues[n]->pos,c->residues[n]);
    }  
    return(0);
  }
}
