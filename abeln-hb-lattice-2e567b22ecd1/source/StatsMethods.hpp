/**
 * \file   StatsMethods.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines inline methods for class Stats
 */

#ifndef _Stats_Methods_H_
#define _Stats_Methods_H_


#include "Stats.hpp"
#include "Lattice.hpp"
#include "Pos.hpp"
#include "AA.hpp"
#include "Chain.hpp"
#include "Native.hpp"
#include "App.hpp"
extern Pos cross[12];

/////////////////////////////////////////////////////////////////////////
inline
void Stats::clean(){
  Esol=0;
  Eext=0;
  Eint=0;
  Etot=0;
  Cint=0;
  Cext=0;
  Ctot=0;
  if(App::nativeStats){
    Nint=0;
    Next=0;
    Ntot=0;
  }

  if(App::hbondStats){
    Hint=0;
    Hext=0;
    Htot=0;
  }

}


inline
Stats::Stats(){
  Esol=0;
  Eext=0;
  Eint=0;
  Etot=0;
  Cint=0;
  Cext=0;
  Ctot=0;
  Nint=0;
  Next=0;
  Ntot=0;
  Hint=0;
  Hext=0;
  Htot=0;

};



inline
Stats Stats::operator+ (const Stats& s) const{
  Stats result;
  result.Cint = (Cint + s.Cint);
  result.Cext = (Cext + s.Cext);
  result.Ctot = (Ctot + s.Ctot);
  result.Eint = (Eint + s.Eint);
  result.Eext = (Eext + s.Eext);
  result.Etot = (Etot + s.Etot);
  result.Esol = (Esol + s.Esol);
  if(App::nativeStats){
    result.Nint = (Nint + s.Nint);
    result.Next = (Next + s.Next);
    result.Ntot = (Ntot + s.Ntot);
  }
  if(App::hbondStats){
    result.Hint = (Hint + s.Hint);
    result.Hext = (Hext + s.Hext);
    result.Htot = (Htot + s.Htot);
  }

  return result;
}


inline
Stats & Stats::operator+= (const Stats& s){
  Stats result;
  Cint += ( s.Cint);
  Cext += ( s.Cext);
  Ctot += ( s.Ctot);
  Eint += ( s.Eint);
  Eext += ( s.Eext);
  Etot += ( s.Etot);
  Esol += ( s.Esol);
  if(App::nativeStats){
    Nint += ( s.Nint);
    Next += ( s.Next);
    Ntot += ( s.Ntot);
  }
  if(App::hbondStats){
    Hint += (s.Hint);
    Hext += (s.Hext);
    Htot += (s.Htot);
  }
  return (*this);
}

inline
Stats & Stats::operator-= (const Stats& s){
  Stats result;
  Cint -= ( s.Cint);
  Cext -= ( s.Cext);
  Ctot -= ( s.Ctot);
  Eint -= ( s.Eint);
  Eext -= ( s.Eext);
  Etot -= ( s.Etot);
  Esol -= ( s.Esol);
  if(App::nativeStats){
    Nint -= ( s.Nint);
    Next -= ( s.Next);
    Ntot -= ( s.Ntot);
  }
  if(App::hbondStats){
    Hint -= (s.Hint);
    Hext -= (s.Hext);
    Htot -= (s.Htot);
  }
  return (*this);
}

inline
Stats & Stats::operator=(const Stats & s){
  Cint =  s.Cint;
  Cext =  s.Cext;
  Ctot =  s.Ctot;
  Eint =  s.Eint;
  Eext =  s.Eext;
  Etot =  s.Etot;
  Esol =  s.Esol; 
  if(App::nativeStats){
    Nint =  s.Nint;
    Next =  s.Next;
    Ntot =  s.Ntot ;
  }

  if(App::hbondStats){
    Hint = (s.Hint);
    Hext = (s.Hext);
    Htot = (s.Htot);
  }

  return (*this);
}

inline
const bool  Stats::operator != (const Stats& s) const{
  bool ans = (Cint !=  s.Cint ||
	      Cext !=  s.Cext ||
	      Ctot !=  s.Ctot ||
	      Eint !=  s.Eint ||
	      Eext !=  s.Eext ||
	      Esol !=  s.Esol ||
	      Etot !=  s.Etot );
  if(App::nativeStats){
    ans = (ans ||
	   Nint !=  s.Nint ||
	   Next !=  s.Next ||
	   Ntot !=  s.Ntot );
  }

  if(App::hbondStats){
    ans = (ans ||
	   Hint !=  s.Hint ||
	   Hext !=  s.Hext ||
	   Htot !=  s.Htot );
  }

  return ans;
}

inline
Stats  Stats::delta(const Stats & Snew, const  Stats & Sold)const{  
  Stats out; 
  out.Eext = Eext + Snew.Eext - Sold.Eext;
  out.Eint = Eint + Snew.Eint - Sold.Eint;
  out.Etot = Etot + Snew.Etot - Sold.Etot;
  out.Esol = Esol + Snew.Esol - Sold.Esol;

  out.Cext = Cext + Snew.Cext - Sold.Cext;
  out.Cint = Cint + Snew.Cint - Sold.Cint;
  out.Ctot = Ctot + Snew.Ctot - Sold.Ctot;

  if(App::nativeStats){	 
    out.Next = Next + Snew.Next - Sold.Next;
    out.Nint = Nint + Snew.Nint - Sold.Nint;
    out.Ntot = Ntot + Snew.Ntot - Sold.Ntot;
  }
  if(App::hbondStats){	 
    out.Hext = Hext + Snew.Hext - Sold.Hext;
    out.Hint = Hint + Snew.Hint - Sold.Hint;
    out.Htot = Htot + Snew.Htot - Sold.Htot;
  }


  return out;
}




inline 
void Stats::solventStats(Pos & pos,Lattice * l){
  if(AA::WATER < 0)return;
  int Es=0;
  for(int k=0;k<6;k++){
    Pos posNB =local[k]+ pos;
    posNB.periodicBoundary();
    Residue * resNB = l->getResidue(posNB);
    if(resNB!=NULL){
      Pos dirSpin = posNB + resNB->spin;
      dirSpin.periodicBoundary();
      // amino acid specific solvent interaction
      if(dirSpin==pos){
	Es +=   l->aaInt->getInteraction(resNB->aa, AA::WATER);
      }else{
	
	/*Pos antiSpin = posNB - resNB->spin;
	antiSpin.periodicBoundary();
	if(antiSpin != pos){
	  //backbone solvent interaction
	  if(resNB->fwd == resNB->spin){
	    Pos antiBkwd = posNB - resNB->fwd;
	    antiBkwd.periodicBoundary();
	    if(antiBkwd != pos)Es += App::BBSolvE;
	  }else if(resNB->fwd== resNB->spin){
	    Pos antiFwd = posNB - resNB->bkwd;
            antiFwd.periodicBoundary();
            if(antiFwd!= pos)Es += App::BBSolvE;
	  }else{
	    Es += App::BBSolvE; 
	  }
	  }*/ 

	// calculate hbond direction (normal) vector
	Pos hbond;
        hbond.periodicSubtraction(pos,posNB);

	// test for orthogonality with spin        
	if(Pos::orthogonal(hbond,resNB->spin)){

	  // test if spin orthogonal on backbone
	  if(Pos::orthogonal(resNB->fwd,resNB->spin) && 
	     Pos::orthogonal(resNB->bkwd,resNB->spin)){
	    // spin orthogonal to backbone
	    // orthoganilty of spin and hbond is sufficient 
	    Es += App::BBSolvE; 
	  }else{
	    // spin in line with backbone
	    // also need hbond orthogonality with backbone
            if(Pos::orthogonal(hbond,resNB->fwd) &&
               Pos::orthogonal(hbond,resNB->bkwd)){
              Es += App::BBSolvE;
            }
	  }
	}

	// alternative implementation
	/*Pos hbond; 
	hbond.periodicSubtraction(pos,posNB);
	if(Pos::orthogonal(hbond,resNB->spin) &&
	   Pos::orthogonal(hbond,resNB->fwd) &&
	   Pos::orthogonal(hbond,resNB->bkwd)){
	   Es += App::BBSolvE;
	   }*/
      }
      /*Pos nbOppSpin; nbOppSpin.periodicSubtraction(posNB, resNB->spin);
	nbOppSpin.periodicBoundary();
	if(AA::WATER > 0 && nbOppSpin!= pos){
	Es +=   l->aaInt->getInteraction(resNB->aa, AA::WATER);
	}*/
    }
  }
  Esol += Es;
  Etot += Es;
}


inline
void Stats::localStats(Residue * r, Lattice *l){
  localStats(r,r->pos,r->spin,r->state,r->aa,l);
}

inline
void Stats::localStats(Residue * r, Pos p,Lattice *l){
  localStats(r,p,r->spin,r->state,r->aa,l);
}


inline
void Stats::localStats(Residue * res,const Pos & pos, const Pos & spin,State state, Lattice * l){
  localStats(res,pos,spin,state,res->aa,l);
}

inline
void Stats::localStats(Residue * res,const Pos & pos, const Pos & spin, State state,int resAA, Lattice * l){
  localStats(res,pos,spin,res->fwd,res->bkwd,state,resAA,l);
}


inline
void Stats::localStats(Residue * res,const Pos & pos, const Pos & spin, const Pos & fwd, const Pos & bkwd, State state,int resAA, Lattice * l){ 
  // for Int(AA) and h-bond energy
  for(int k=0;k<6;k++){
    Pos posNB =local[k]+ pos;
    posNB.periodicBoundary();
    Residue * resNB = l->getResidue(posNB);
    if(resNB!=NULL && resNB != l->AIR){
      bool sameChain = resNB->chainNum==res->chainNum;
      if(abs(resNB->n - res->n)!=1 || ! sameChain){
	// Neighbouring sites on lattice, excluding neighbours in chain
	// only count non-masked residues 
	if(!res->masked &&  !resNB->masked){ 
	  if(sameChain){Cint++;}else{Cext++;}
	}	  
	if(App::nativeStats){
	  if(l->native->isNativeContact(res->chainNum,resNB->chainNum,res->n,resNB->n)){
	  // only count non-masked residues 
	  if(!res->masked &&  !resNB->masked){
	    if(sameChain){Nint++;}else{Next++;}
	  }
	}
      }

	Pos fromRes = spin+pos;
	fromRes.periodicBoundary();
	Pos fromNB =resNB->pos + resNB->spin;
	fromNB.periodicBoundary();
	bool spinsAligned= resNB->spin==spin &&  fromNB != pos && fromRes != resNB->pos;
	//hbonds:
	if( state == inStrand && resNB->state == inStrand && spinsAligned){
	  //inStrand && alligned
	  if(sameChain){
	    Eint += App::HBondE;
	  }else{
	    Eext += App::HBondE;
	  }
	  if(App::hbondStats){
	  // only count non-masked residues 
	  if(!res->masked &&  !resNB->masked){
	    if(sameChain){Hint++;}else{Hext++;}
	  }
	}      
      }      
	//AA interactions, depending on Spin
	bool spinsFacing= (resNB->spin == (Pos(0,0,0)- spin) && fromNB== pos);
	if(spinsAligned || spinsFacing){
	  if(sameChain){
	    Eint +=  l->aaInt->getInteraction(resAA, resNB->aa);
	  }else{
	    Eext +=  l->aaInt->getInteraction(resAA, resNB->aa);
	  }
	}
      }else{ // neighbours in Chain 
	//check steric hindrance
	if(sameChain && spin == resNB->spin){
	  Eint += App::StericE;
	}
      }
    }else if(resNB != l->AIR){ // water contact
    // Pos nbOppSpin; nbOppSpin.periodicSubtraction(pos, spin);
    // nbOppSpin.periodicBoundary();
    // if(AA::WATER > 0 && nbOppSpin!= posNB){
    Pos dirSpin = pos + spin;
    dirSpin.periodicBoundary();
    if(AA::WATER >= 0){
      if( posNB == dirSpin){
	Esol +=   l->aaInt->getInteraction(resAA, AA::WATER);
      }else{
	// calculate h-bond interaction with solvent
	
	// calculate hbond direction
	Pos hbond;                                                                          
        hbond.periodicSubtraction(pos,posNB);  
	
	// test for orthogonality with spin                                                        
        if(Pos::orthogonal(hbond,spin)){

          // test if spin orthogonal with backbone       
	  if(Pos::orthogonal(fwd,spin) &&
             Pos::orthogonal(bkwd,spin)){
            // spin orthogonal to backbone                                          
            // orthoganilty of spin and hbond is sufficient       
            Esol += App::BBSolvE;
          }else{
            // spin in line with backbone              
            // also need hbond orthogonality with backbone                                       
            if(Pos::orthogonal(hbond,fwd) &&
               Pos::orthogonal(hbond,bkwd)){
              Esol += App::BBSolvE;
            }
          }
	}

	// alternative implementation - wrong, can have 3 hbonds on corner
	/*	Pos hbond;
        hbond.periodicSubtraction(pos,posNB);
        if(Pos::orthogonal(hbond,spin)&&
	   Pos::orthogonal(hbond,fwd) && 
	   Pos::orthogonal(hbond,bkwd)){
	  Esol += App::BBSolvE;
	  }*/
      }
    }
    }
  }

#ifdef CROSSNBS
Pos * crossNbs = getCrossNbsSpinDir(spin);
for(int k=0;k<4;k++){
  Pos posNB =crossNbs[k]+ pos;
  posNB.periodicBoundary();
  Residue * resNB = l->getResidue(posNB);
  if(resNB!=NULL && resNB != l->AIR){
    bool sameChain = resNB->chainNum==res->chainNum;
     
    //only facing spin interactions...
    bool spinsOpposite= (resNB->spin == (Pos(0,0,0)- spin));
    if(spinsOpposite){
      if(sameChain){
	Eint +=  l->aaInt->getInteraction(resAA, resNB->aa)/ CROSS_FACTOR;
	Cint++;
      }else{
	Eext += l->aaInt->getInteraction(resAA, resNB->aa)/CROSS_FACTOR;
	Cext++;
      }
    }
  }else if( resNB != l->AIR){ // water contact
    /*Pos nbOppSpin; nbOppSpin.periodicSubtraction(pos, spin);
      nbOppSpin.periodicBoundary();
      if(AA::WATER > 0 && nbOppSpin!= posNB){
      Esol +=   l->aaInt->getInteraction(res->aa, AA::WATER);
      }*/
    Pos dirSpin = pos + spin;
    dirSpin.periodicBoundary();
    if(AA::WATER >= 0 && posNB == dirSpin){
      Esol +=   l->aaInt->getInteraction(resAA, AA::WATER);
    }

  }
 }
//delete crossNbs;
#endif
Etot=Eext+Eint+Esol;
Ctot=Cint+Cext;
if(App::nativeStats){
  Ntot=Nint+Next;
 }
if(App::hbondStats){
  Htot=Hint+Hext;
 }
}

inline
void Stats::localStatsExclude(Residue * res,const Pos & pos, const Pos & spin,State state,
			      Lattice * l,int start,int end){
  
  // for Int(AA) and h-bond energy
  for(int k=0;k<6;k++){
    Pos posNB =local[k]+ pos;
    posNB.periodicBoundary();
    Residue * resNB = l->getResidue(posNB);
    if(resNB!=NULL && resNB != l->AIR){
      bool sameChain = resNB->chainNum==res->chainNum;
      if(!sameChain ||( resNB->n <start || resNB->n > end)){ // exluded residues
	if(abs(resNB->n - res->n)!=1 || ! sameChain){
	  // Neighbouring sites on lattice, excluding neighbours in chain
	  
	  // only count non-masked residues 
	  if(!res->masked &&  !resNB->masked){
	    if(sameChain){Cint++;}else{Cext++;}
	  }	  

	  if(App::nativeStats){
	    if(l->native->isNativeContact(res->chainNum,resNB->chainNum,res->n,resNB->n)){
	     // only count non-masked residues 
	      if(!res->masked &&  !resNB->masked){
		if(sameChain){Nint++;}else{Next++;}
	      }
	    }
	  }
	  Pos fromRes =pos + spin;
	  fromRes.periodicBoundary();
	  Pos fromNB =resNB->pos + resNB->spin;
	  fromNB.periodicBoundary();
	  bool spinsAligned= resNB->spin== spin &&  fromNB != pos && fromRes != resNB->pos;
	  //hbonds:
	  if( state == inStrand && resNB->state == inStrand && spinsAligned){
	    //inStrand && alligned
	    if(sameChain){
	      Eint += App::HBondE;
	    }else{
	      Eext += App::HBondE;
	    }
	    if(App::hbondStats){
	      // only count non-masked residues 
	      if(!res->masked &&  !resNB->masked){
		if(sameChain){Hint++;}else{Hext++;}
	      }
	    }
	  }      
	  //AA interactions, depending on Spin
	  bool spinsFacing= (resNB->spin == (Pos(0,0,0)- spin)) && fromNB== pos ;
	  if(spinsAligned || spinsFacing){
	    if(sameChain){
	      Eint +=  l->aaInt->getInteraction(res->aa, resNB->aa);
	    }else{
	      Eext +=  l->aaInt->getInteraction(res->aa, resNB->aa);
	    }
	  }
	}else{ // neighbours in Chain 
	  //check steric hindrance
	  if(sameChain && spin == resNB->spin){
	    Eint += App::StericE;
	  }
	}
      }
    }else if( resNB != l->AIR){ // water contact
      /*Pos nbOppSpin; nbOppSpin.periodicSubtraction(pos, spin);
	nbOppSpin.periodicBoundary();
	if(AA::WATER > 0 && nbOppSpin!= posNB){
	Esol +=   l->aaInt->getInteraction(res->aa, AA::WATER);
	}*/
      Pos dirSpin = pos + spin;
      dirSpin.periodicBoundary();
      if(AA::WATER >= 0){
	if( posNB == dirSpin){
	  Esol +=   l->aaInt->getInteraction(res->aa, AA::WATER);
	}else{
	  Pos antiSpin = pos - spin;
	  antiSpin.periodicBoundary();
	  if(antiSpin != pos){
	    //backbone solvent interaction                                                                         	    Esol +=App::BBSolvE;
	  }
	}
      }
    }
  }
#ifdef CROSSNBS
  Pos * crossNbs = getCrossNbsSpinDir(spin);
  for(int k=0;k<4;k++){
    Pos posNB =crossNbs[k] + pos;
    posNB.periodicBoundary();
    Residue * resNB = l->getResidue(posNB);
    if(resNB!=NULL && resNB != l->AIR ){
      bool sameChain = resNB->chainNum==res->chainNum;
      if(!sameChain ||( resNB->n <start || resNB->n > end)){
	//only facing spin interactions...
	bool spinsOpposite= (resNB->spin == (Pos(0,0,0)- spin));
	if(spinsOpposite){
	  if(sameChain){
	    Eint += l->aaInt->getInteraction(res->aa, resNB->aa)/CROSS_FACTOR;
	    Cint++;
	  }else{
	    Eext += l->aaInt->getInteraction(res->aa, resNB->aa)/CROSS_FACTOR;
	    Cext++;
	  }
	}
      }
    }else if( resNB != l->AIR){ // water contact
      /*Pos nbOppSpin; nbOppSpin.periodicSubtraction(pos, spin);
	nbOppSpin.periodicBoundary();
	if(AA::WATER > 0 && nbOppSpin!= posNB){
	Esol +=   l->aaInt->getInteraction(res->aa, AA::WATER);
	}*/
      Pos dirSpin = pos + spin;
      dirSpin.periodicBoundary();
      if(AA::WATER >= 0 && posNB == dirSpin){
        Esol +=   l->aaInt->getInteraction(res->aa, AA::WATER);
      }
    }
  }
  // delete crossNbs;
#endif
  Etot=Eext+Eint+Esol;
  Ctot=Cint+Cext;
  if(App::nativeStats){
    Ntot=Nint+Next;
  }
  if(App::hbondStats){
    Htot=Hint+Hext;
  }
}




inline
int Stats::getDeltaE(const Stats & old)const{
  return Etot -old.Etot;
}




//static
#ifdef CROSSNBS
inline
Pos * Stats::getCrossNbsSpinDir(Pos spin){
  //cout<<"spin "<<spin.toString()<<endl;
  //cout<<loopUpCrossNbs[spin.x+1][spin.y+1][spin.z+1][0].toString()<<endl;
  return loopUpCrossNbs[spin.x+1][spin.y+1][spin.z+1];
}
#endif



/*inline
//static
Pos * Stats::getCrossNbsSpinDir(Pos spin){
Pos *ans= new Pos[4];
int indx=0;
for(int i=0;i<12;i++){
if((spin.x !=0 && cross[i].x==spin.x) || 
(spin.y !=0 && cross[i].y==spin.y) ||
(spin.z !=0 && cross[i].z==spin.z)){
ans[indx]=cross[i];
indx++;
}
}
if(indx != 4){cout <<"too few neighbours"<<indx<<endl;exit(1);}
return ans;
}*/



#endif


