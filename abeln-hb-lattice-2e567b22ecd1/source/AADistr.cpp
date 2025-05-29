#include "AADistr.hpp"
#include "AA.hpp"
#include "Lattice.hpp"
#include "Residue.hpp"
#include "Chain.hpp"

#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>


AADistr::AADistr(){

}

void AADistr::init(string filename, AA * aai){
  aaInt =aai;
  numaa = AA::NUMAA;
  sumFreq=0;
  //// read file
  ifstream aaFile (filename.c_str());
  if (!aaFile.is_open()){
    cout << "could not open file: "<<filename.c_str()<<  endl;
    exit(1);
  }
  int i=0;


  char line[500];
  while(aaFile.getline(line,200)){
    
    
    
    //// split line into fields
    char aa_c[5];
    int freq;
    sscanf(line,"%3s%d",aa_c,&freq);
    string aa = string(aa_c);
    
    cout<< aa<<endl;
    //// check if aa order is same as in neteraction matrix
    if(aa != "HOH"){
      if( aaInt->stringtoAA(aa)!=i){
	cout<< "amino acid order in file "<<filename<<" does not match interaction matrix";
	exit(1);
      }
      freqPDB[i]=freq;
      sumFreq+=freq;
    }
    i++;
  }  
}

void AADistr::setLength(int length){
  seqLength=length;
  for(int i=0;i<numaa;i++){
    if(i != AA::WATER){ 
      freqExp[i]=freqPDB[i]*seqLength/sumFreq;
    }else{
      freqExp[i]=0;
    }
  }

}


void AADistr::setSequence(Lattice *l){
  int length=0;
  //// reset frequency table
  for(int i=0;i<AA::NUMAA;i++){
    freqSeq[i]=0;
  }
  //// 
  for(int nc=0;nc<l->nChains;nc++){
    Chain * ch = l->chains[nc];
    for (int n=0;n<ch->N;n++){
      Residue * res =ch->residues[n];
      freqSeq[res->aa]++;
      length++;
    }
  }
  setLength(length);
  systemDistance = getDistance();
}




double AADistr::getDistance(){
  double dist=0;
  for(int i=0;i<numaa;i++){
    if(i != AA::WATER){ 
      dist += (freqSeq[i] - freqExp[i])*(freqSeq[i] - freqExp[i]);
    }
  }
  return dist;
}


bool AADistr::check(){
  double calcD=getDistance();
  if(calcD != systemDistance ){
    cout <<"ERROR in design process accummulative distance:"<<systemDistance;
    cout <<" does not match calc distance:" <<calcD<<endl;
    return false;
  }
  return true;
}

double myround(double x,int decimal_places){
  double tmp = pow((double)10.0,(double)decimal_places)*x;
  return floor(tmp + 0.5)/  pow((double)10,(double)decimal_places);
}

string AADistr::toString(){
  std::stringstream outs;
  outs<<"system distance " <<systemDistance<<endl;
  outs<<"SEQUENCE"<<endl;
  double l=0;
  for(int aa=0;aa<AA::NUMAA;aa++){
    outs <<AA::int2aa[aa]<<" "<<freqSeq[aa]<<", ";
    l+= freqSeq[aa];
  }
  outs<< " length " <<l<< endl;
  outs<<"PDB"<<endl;
  double checksum =0;
  for(int aa=0;aa<AA::NUMAA;aa++){
    outs <<AA::int2aa[aa]<<" "<<myround(freqExp[aa],1)<<", ";
    checksum += freqExp[aa];
  }
  outs<< "SUM: "<<checksum <<endl;
  return outs.str();
} 

