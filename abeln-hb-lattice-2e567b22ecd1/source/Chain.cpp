#include "Chain.hpp"
#include "Lattice.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "EnergyMap.hpp"

#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;


Chain::Chain(Lattice *l1,int chnum){
  l=l1;
  chainNum=chnum;
  frozen =false;
  locked=false;

}


Chain::~Chain(){
  for(int i=0;i<N;i++){
    delete residues[i];
  }

}

void Chain::setChainNum(int n){
  chainNum=n;
}


void Chain::setRandomSpins(){
  setFwdBkwdSpins();
  Pos dir;
  for(int n=0;n<N;n++){
    do{
      dir = local[(int) floor(6 * drand48())];
    }while(dir == residues[n]->bkwd || dir == residues[n]->fwd );
    residues[n]->spin= dir;
  }  
}

void Chain::setFwdBkwdSpins(){

  residues[0]->bkwd=Pos(0,0,0);

  residues[0]->fwd.periodicSubtraction(residues[1]->pos,residues[0]->pos);
 
  residues[N-1]->fwd=Pos(0,0,0);
  residues[N-1]->bkwd.periodicSubtraction(residues[N-2]->pos, residues[N-1]->pos);
 
  for(int n=1;n<N-1;n++){
    residues[n]->bkwd.periodicSubtraction(residues[n-1]->pos,residues[n]->pos);
    residues[n]->fwd.periodicSubtraction(residues[n+1]->pos,residues[n]->pos);
  }
}

void Chain::setFwdBkwd(Residue * res){
  int n= res->n;
  
  if(n!=0){
    res->bkwd.periodicSubtraction(residues[n-1]->pos, residues[n]->pos);
  }else{
    res->bkwd = Pos(0,0,0);
  }
  
  if(n != (N-1)){
    res->fwd.periodicSubtraction(residues[n+1]->pos,residues[n]->pos);
  }else{
    res->fwd = Pos(0,0,0);
  }
  
}


bool Chain::checkFwdBkwd(Residue * res){
  return (res->spin != res->bkwd && res->spin != res->fwd);
}

bool Chain::checkAllFwdBkwd(){
  for(int i=0;i<N;i++){
    if(!checkFwdBkwd(residues[i])){
      cout<<"SPIN FWD BKWD VIOLATION Residue:"<<i<< "chain: "<< chainNum<< endl;
      cout<<"fwd : "<<residues[i]->fwd.toString() <<endl;
      cout<<"bkwd: "<<residues[i]->bkwd.toString()   <<endl;
      cout<<"spin: "<<residues[i]->spin.toString()   <<endl;
      l->writePDBMultiChain("spinViolation.pdb");
      exit(1);
    }
  }
  return true;
}


void Chain::setStateResidues(){
  for(int n=0;n<N;n++){
    residues[n]->state=inCoil;
  }
}




void Chain::newRandomSpin(Residue * res){
  do{
    res->spin = local[(int) floor(6 * drand48())];
  }while(res->spin == res->bkwd || res->spin == res->fwd );
}


Pos Chain::newSpinPosCornerFlip(Pos oldSpin,Pos n1, Pos n2){
  //cout<<"spinCornerFlipStarted"<<endl;
  //cout <<"n1 "<<n1.toString()<<endl;
  //cout <<"n2 "<<n2.toString()<<endl;
  //cout <<"oldSpin "<<oldSpin.toString()<<endl;
  Pos diff;
  diff.periodicSubtraction(n1,n2);
  //cout <<"diff "<<diff.toString()<<endl;
  int pp=-1,o1=9,o2=-1;
  for(int dir=0;dir<3;dir++){
    if(diff.xyz[dir]==0){pp=dir;}
    else if(o1==9){o1=dir;}
    else{o2=dir;}
  }
  // spin othogonal to plane of move
  Pos result;//initialized to zero
  if(oldSpin.xyz[pp]!=0){ return (result-oldSpin );}
  // check direction of move
  if(diff.xyz[o1]==diff.xyz[o2]){
    result.xyz[o1]=oldSpin.xyz[o2];
    result.xyz[o2]=oldSpin.xyz[o1];
  }else{
    result.xyz[o1]=  - oldSpin.xyz[o2];
    result.xyz[o2]=  - oldSpin.xyz[o1];
  }
  return result;
}

Pos Chain::newSpinEndMove(Pos oldPos, Pos newPos, Pos posTail, Pos oldSpin){
  //cout<<"spinEndMoveStarted"<<endl;
  //cout <<"oldPos "<<oldPos.toString()<<endl;
  //cout <<"newPos "<<newPos.toString()<<endl;
  // cout <<"posTail "<<posTail.toString()<<endl;
  //cout <<"oldSpin "<<oldSpin.toString()<<endl;
  Pos diff;
  diff.periodicSubtraction(oldPos,newPos);
  //cout <<"diff "<<diff.toString()<<endl;
  int pp=-1,o1=9,o2=-1;
  int cntpp=0;
  for(int dir=0;dir<3;dir++){
    if(diff.xyz[dir]==0){pp=dir;cntpp++;}
    else if(o1==9){o1=dir;}
    else{o2=dir;}
  }


  if(cntpp>1){
    //180 degree flip
    Pos tail;
    tail.periodicSubtraction(oldPos,posTail);
    int dirtail=0;
    for(int dir=0;dir<3;dir++){
      if(tail.xyz[dir]!=0){dirtail=dir;}
    }
    if(oldSpin.xyz[o1]==0 && oldSpin.xyz[dirtail]==0){
      //orthogonal to moving plane
      return oldSpin;
    }else{
      return (Pos(0,0,0) - oldSpin);
    }
  }

 
  //cout<<pp<<" "<<o1<<" "<<o2<<endl;   

  // spin orthogonal  to plane of move
  if (oldSpin.xyz[pp]!=0){
    return oldSpin;
  }
  Pos tail;
  tail.periodicSubtraction(oldPos,posTail);
  if(tail == oldSpin){
    Pos newTail;
    newTail.periodicSubtraction(newPos,posTail);
    return newTail;
  }
  

  Pos result;
  if (diff.xyz[o1]!=diff.xyz[o2]){
   
    result.xyz[o1]=oldSpin.xyz[o2];
    result.xyz[o2]=oldSpin.xyz[o1];

  }else{
  
    result.xyz[o1]= - oldSpin.xyz[o2];
    result.xyz[o2]= - oldSpin.xyz[o1];
  }
  //cout<<"spinEndMoveDone"<<endl;
  return result;
}


bool Chain::strandPossible(Residue * res){
  int n = res->n;
  if((n==0) || (n==N-1)){
    return false;
  }else{
    if(res->fwd != (Pos(0,0,0)- res->bkwd)){
      return false;
    }
  }
  if(n>0){
    Residue * prev=residues[n-1];
    if(prev->state==inStrand){
      if(prev->spin !=   (Pos(0,0,0)- res->spin)){
	return false;
      }
    }
  }
  if(n<(N-1)){
    Residue * nxt=residues[n+1];
    if(nxt->state==inStrand){
      if(nxt->spin !=   (Pos(0,0,0)- res->spin)){
	return false;
      }
    }
  }
  return true;
}



int Chain::getCext(){
  int Cext=0;
  for(int n=0;n<N;n++){
    Residue * res= residues[n];
    //if(res->aa != 13){
    {
      for(int k=0;k<6;k++){
	Pos posNB = local[k] + res->pos;
        posNB.periodicBoundary();
	Residue * neighbour =l->getResidue(posNB); 
	if(neighbour!=NULL && neighbour->chainNum != res->chainNum){
	  //if((neighbour->aa != 13)) Cext++;
	  Cext++;
	}
      }
#ifdef CROSSNBS
      Pos * crossNbs = Stats::getCrossNbsSpinDir(res->spin);
      for(int k=0;k<4;k++){
	Pos posNB =crossNbs[k]+ res->pos;
	posNB.periodicBoundary();
	Residue * resNB = l->getResidue(posNB);
	if(resNB!=NULL){
	  bool sameChain = resNB->chainNum==res->chainNum;
	  if(!sameChain){
	    bool spinsOpposite = resNB->spin == (Pos(0,0,0)- res->spin);
	    if(spinsOpposite){
	      Cext++;
	    }
	  }
	}
      }
      //delete crossNbs;
#endif
    }
  }
  return Cext;
}


bool Chain::hasCext(){
  for(int n=0;n<N;n++){
    Residue * res= residues[n];
    //if(res->aa != 13){
    
    for(int k=0;k<6;k++){
      Pos posNB = local[k] + res->pos;
      posNB.periodicBoundary();
      Residue * neighbour =l->getResidue(posNB); 
      if(neighbour!=NULL && neighbour->chainNum != res->chainNum){
	//if((neighbour->aa != 13)) 
	return true;
      }
    }
    
#ifdef CROSSNBS
    Pos * crossNbs = Stats::getCrossNbsSpinDir(res->spin);
    for(int k=0;k<4;k++){
      Pos posNB =crossNbs[k]+ res->pos;
      posNB.periodicBoundary();
      Residue * resNB = l->getResidue(posNB);
      if(resNB!=NULL){
	bool sameChain = resNB->chainNum==res->chainNum;
	if(!sameChain){
	  bool spinsOpposite= (resNB->spin == (Pos(0,0,0)- res->spin));
	  if(spinsOpposite){
	    return true;
	  }
	}
      }
    }
    //delete crossNbs;
#endif
  }
  return false;
}





Pos Chain::rotationPosition(Pos posOldi, int rotationDir,Pos posPivot){
  Pos output;
  int prim = rotationDir % 3;
  int sign =1;
  if(rotationDir >2){
    sign=-1;
  }
  Pos posOld =posOldi - posPivot ;
  switch(prim){
  case 0: // around x-axis
    output.x = posOld.x+ posPivot.x;
    output.y = sign*posOld.z + posPivot.y ;
    output.z= -1*sign*posOld.y + posPivot.z;
    break;
  case 1:
    output.x = -1*sign*posOld.z + posPivot.x;
    output.y = posOld.y + posPivot.y;
    output.z = sign*posOld.x + posPivot.z;
    break;
  case 2:
    output.x = sign*posOld.y + posPivot.x;
    output.y = -1*sign*posOld.x + posPivot.y ;
    output.z = posOld.z + posPivot.z;
    break;
  default:
    cout << "rotationCoord: wrong direction: "<< prim<<endl;
    exit(1);
  }
  output.periodicBoundary();
  return output  ;
}



Pos Chain::rotationSpin(Pos spinOld, int rotationDir){
  Pos output;
  int prim = rotationDir % 3;
  int sign =1;
  if(rotationDir >2){
    sign=-1;
  }
  switch(prim){
  case 0: // around x-axis
    output.x = spinOld.x ;
    output.y = sign*spinOld.z;
    output.z= -1*sign*spinOld.y ;
    break;
  case 1:
    output.x = -1*sign*spinOld.z ;
    output.y = spinOld.y ;
    output.z = sign*spinOld.x ;
    break;
  case 2:
    output.x = sign*spinOld.y ;
    output.y = -1*sign*spinOld.x  ;
    output.z = spinOld.z;
    break;
  default:
    cout << "rotationCoord: wrong direction: "<< prim<<endl;
    exit(1);
  }
  return output  ;
}




int Chain::rotationCoord(int x,int y,int z, int dir, int rotationDir,int pivot){
  // OPT output 3D vector
  int prim = rotationDir % 3;
  int sign =1;
  if(rotationDir >2){
    sign=-1;
  }
  x -= residues[pivot]->pos.x;
  y -= residues[pivot]->pos.y;
  z -= residues[pivot]->pos.z;

  switch(prim){
  case 0: // around x-axis
    if(dir==XDIR) return x+ residues[pivot]->pos.x;
    else if(dir==YDIR) return ((sign*z + residues[pivot]->pos.y)  +LY)%LY;
    else return ((-1*sign*y + residues[pivot]->pos.z)   +LZ)%LZ;
    //break;
  case 1: // around y-axis
    if(dir==XDIR) return ((-1*sign*z + residues[pivot]->pos.x)+ LX)%LX;
    else if(dir==YDIR) return y + residues[pivot]->pos.y;
    else return ((sign*x + residues[pivot]->pos.z) + LZ)%LZ;
    //break;
  case 2: // around z-axis
    if(dir==XDIR) return ((sign*y + residues[pivot]->pos.x) +LX)%LX;
    else if(dir==YDIR) return ((-1*sign*x + residues[pivot]->pos.y) +LY)%LY;
    else return z + residues[pivot]->pos.z;
    //break;
  default:
    cout << "rotationCoord: wrong direction: "<< prim<<endl;
    exit(1);
   }
  return(-1);
}






