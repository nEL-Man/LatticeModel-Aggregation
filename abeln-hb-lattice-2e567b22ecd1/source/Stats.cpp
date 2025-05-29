#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "Lattice.hpp"
#include "Native.hpp"
#include "AA.hpp"


#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

extern Pos cross[12];

#ifdef CROSSNBS
Pos Stats::loopUpCrossNbs[3][3][3][4];
#endif

int getInteraction(int,int);

#ifdef CROSSNBS
void Stats::setLoopUpCrossNbs(){
  int cnt[3][3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	cnt[i][j][k]=0;
      }
    }
  }
  for(int k=0;k<6;k++){
    Pos spin = local[k];
    for(int i=0;i<12;i++){
      if((spin.x !=0 && cross[i].x==spin.x) || 
	 (spin.y !=0 && cross[i].y==spin.y) ||
	 (spin.z !=0 && cross[i].z==spin.z)){
	int indx = cnt[spin.x+1][spin.y+1][spin.z+1];
	loopUpCrossNbs[spin.x+1][spin.y+1][spin.z+1][indx]=cross[i];
	cnt[spin.x+1][spin.y+1][spin.z+1]++;
      }    
    }
  }
  return;
}
#endif

void Stats::getLatticeStats(Lattice * l){
  //cout<< "getting Lattice stats"<<endl;
  Eext=0;
  Eint=0;
  Cint=0;
  Cext=0;
  Ctot=0;
  Etot=0;
  Esol=0;	 
  Nint=0;
  Next=0;
  Ntot=0;	 
  Hint=0;
  Hext=0;
  Htot=0;


  for(int cn=0;cn< l->nChains ; cn++){
    Chain * c=l->chains[cn];
    for(int n=0;n< c->N;n++){
      // cout<<n<<endl;
      Residue * res= c->residues[n];
      Stats tmp;
      tmp.localStats(res,res->pos,res->spin,res->state,l);
      (*this) += tmp;
    }
  }
  Eext /=2 ;
  Eint /=2;
  Cint /=2;
  Cext /=2;

  Etot=Eext+Eint+Esol;
  Ctot=Cint+Cext;

  Nint /=2;
  Next /=2;
  Ntot=Nint+Next;
  Hint /=2;
  Hext /=2;
  Htot=Hint+Hext;
}





void Stats::get_Eint_Cint(Chain*  c, Lattice *l){
  Eint =0;
  Cint =0;
  Esol=0;
  for(int n=0;n< c->N;n++){
      Residue * res= c->residues[n];
      Stats tmp;
      tmp.localStats(res,res->pos,res->spin,res->state,l);
      (*this) += tmp;
  }

  Eint /=2;
  Cint /=2;
  Etot = Eint+Esol;
  Ctot = Cint;
  // set all others to zero
  Eext=0;
  Cext=0;

  Nint/=2;  
  Next=0;
  Ntot=Nint;


  Hint/=2;  
  Hext=0;
  Htot=Hint; 
}




void Stats::printCout(){
  cout<<"Eint "<< Eint;
  cout<<" Eext "<< Eext;
  cout<<" Esol "<< Esol;
  cout<<" Etot "<< Etot<<endl;

  cout<<"Cint "<<Cint;
  cout<<" Cext "<<Cext;
  cout<<" Ctot "<<Ctot<<endl;
  

  if(App::nativeStats){	 
    cout<<"Nint "<<Nint;
    cout<<" Next "<<Next;
    cout<<" Ntot "<<Ntot<<endl;
  }
  
  if(App::hbondStats){ 
    cout<<"Hint "<<Hint;
    cout<<" Hext "<<Hext;
    cout<<" Htot "<<Htot<<endl;
  }
}

string Stats::print2string(){
  std::stringstream myStream;
  myStream<<"Eint "<< Eint;
  myStream<<" Eext "<< Eext;
  myStream<<" Esol "<< Esol;
  myStream<<" Etot "<< Etot<<endl;

  myStream<<"Cint "<<Cint;
  myStream<<" Cext "<<Cext;
  myStream<<" Ctot "<<Ctot<<endl;
  

  if(App::nativeStats){	  
    myStream<<"Nint "<<Nint;
    myStream<<" Next "<<Next;
    myStream<<" Ntot "<<Ntot<<endl;
  }

  if(App::hbondStats){   
    myStream<<"Hint "<<Hint;
    myStream<<" Hext "<<Hext;
    myStream<<" Htot "<<Htot<<endl;
  }

  return myStream.str();
}




