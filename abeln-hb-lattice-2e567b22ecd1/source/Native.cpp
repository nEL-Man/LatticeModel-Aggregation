#include "Native.hpp"







Native::Native(Lattice *l){
  cout<<"creating new native"<<endl;
  contacts = new  bool  *** [l->nChains];
  for(int i=0;i<l->nChains;i++){
    contacts[i] = new bool ** [l->nChains];
    int length_i =l->chains[i]->N;
    for(int j=0;j<l->nChains;j++){
      int length_j =l->chains[j]->N;
      contacts[i][j] = new  bool* [length_i];
      for(int k=0;k<length_i;k++){
	contacts[i][j][k]=new bool[length_j];
	for(int l=0;l<length_j;l++){
	  contacts[i][j][k][l]=false;
	}
      }
    }
  }
}




void Native::setNative(Lattice * l){
  cout<<"setting native"<<endl;
  //TO DO should check if not out of bounds
  int totalC=0;
  for(int nc=0;nc<l->nChains;nc++){
    for(int n=0;n<l->chains[nc]->N;n++){
      Residue * res = l->chains[nc]->residues[n];
      for(int k=0;k<6;k++){
	Pos posNB =local[k]+ res->pos;
	posNB.periodicBoundary();
	Residue * resNB = l->getResidue(posNB);
	if(resNB!=NULL){
	  if(resNB->chainNum !=  res->chainNum || abs(resNB->n - res->n) !=1){
	    contacts[res->chainNum][resNB->chainNum][res->n][resNB->n]=true;
	    contacts[resNB->chainNum][res->chainNum][resNB->n][res->n]=true;
	    totalC++;
	  }
	}
      }
    }
  }
  totalC =totalC/ 2;
  totCnat=totalC;
  cout<<"TOTAL NATIVE CONTACTS:" <<totCnat<<endl;
  //l->stats.getLatticeStats(l); //should be set after reading in
  nativeEnergy = l->stats.getEtot();
  cout<<"NATIVE INTERNAL ENERGY:" <<nativeEnergy<<endl;
}
 


