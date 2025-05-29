#include "Lattice.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "EnergyMap.hpp"
#include "AA.hpp"
#include "Native.hpp"
#include "GranCan.hpp"
#include "MolBox.hpp"
#include "App.hpp"

#include <fstream>

Pos local[6]={Pos(1,0,0),
	      Pos(0,1,0),
	      Pos(0,0,1),
	      Pos(-1,0,0),
	      Pos(0,-1,0),
	      Pos(0,0,-1)
};

Pos cross[12]={Pos(1,1,0),
	      Pos(1,-1,0),
	      Pos(-1,1,0),
	      Pos(-1,-1,0),
	      Pos(0,1,1),
	      Pos(0,1,-1),
	      Pos(0,-1,1),
	      Pos(0,-1,-1),
	      Pos(1,0,1),
	      Pos(1,0,-1),
	      Pos(-1,0,1),
	      Pos(-1,0,-1)
};






Lattice::Lattice(){
  //Stats::setLattice(this);
  energyMap = new EnergyMap(this);
  for (int x=0;x<LX;x++){
    for (int y=0;y<LY;y++){
      for (int z=0;z<LZ;z++){
	r[x][y][z]=NULL;
	solv[x][y][z]=water;
      }
    }
  }
  AIR = new Residue;
  AIR->chainNum=-1;
  AIR->n=-1;
}

void Lattice::setSolvent(){
  int xsolv = LX/2;
  for (int y=0;y<LY;y++){
    for (int z=0;z<LZ;z++){
      if( r[xsolv][y][z]!=NULL){
	cout<<"air clashes with chain"<<endl;
	exit(1);
      }
      r[xsolv][y][z]=AIR;
      solv[xsolv][y][z]=air;
    }
  }
}


void Lattice::setAA(string fnaa, int indx_water){
  aaInt=new AA(fnaa,indx_water);
}


void Lattice::readPDBMultiChain(string s){
  //TO DO check if correct format
  cout<<"reading file "<<s<<endl;
  ifstream infile(s.c_str());
  if (!infile.is_open()){
    cout << "can't open file "<< s <<endl;
    exit(0);
  }
  int resnum;
  int cn =0;
  char line[500];
  char type[5],restype[4], atom[4];;
  float xcoord,ycoord,zcoord;
  float occupancy,bfactor;
  // bool lig=false;
  chains[cn]=new Chain(this,cn);
  
  while(!infile.eof()){
    if(cn>= MAX_CHAINS){
      cout<<"Too many chains in "<< s <<endl;
      exit(0);
    }
    infile.getline(line,2000);
    type[0]=' ';
    sscanf(line,"%4s%*9c%2s%*2c%3s%*2c%4d%*4c%9f%9f%9f%4f%4f\n",type,atom,restype,&resnum,&xcoord,&ycoord,&zcoord,&occupancy,&bfactor);
    if ((strcmp(type,"ATOM")==0) &&(strcmp(line,"")!=0)){

      if(resnum>= MAX_RES){
	cout<<"Too many residues in chain "<<cn<< " in file "<< s <<endl;
	exit(0);
      }
      cout<<type<<"*"<<atom<<"*"<<restype<<"*"<<resnum<<"*"<<xcoord<<"*"<<ycoord<<"*"<<zcoord<<"*"<<occupancy<<"*"<<bfactor<<endl;

      if(strcmp(atom,"CA")==0){
	
	Residue * res = new Residue; 
	res->pos.x =(int) (xcoord/3.0);
	res->pos.y = (int) (ycoord/3.0);
	res->pos.z =(int) (zcoord/3.0);
	res->aa = aaInt->stringtoAA(restype);
	res->n =resnum;
	res->chainNum = cn;
	chains[cn]->residues[resnum] = res;
	r[res->pos.x][res->pos.y][res->pos.z]= res;
	if(bfactor==11.0){
	  res->state=inStrand;
	}else if(bfactor==22.0){
	  res->state=inCoil;
	}else{
	  cout<<"state not defined, bfactor: "<<bfactor<<endl;
	  exit(1);
	}
	if(occupancy==1.0){
	  res->masked = false;
	}else if(occupancy==2.0){
	  res->masked = true;
	}else{
	  cout<<"mask state not defined, occupance: "<<occupancy<<endl;
	  exit(1);
	}

      }else if(strcmp(atom,"CB")==0){
	
	Residue * res =chains[cn]->residues[resnum];
        //TO DO periodic boundaries ?
	res->spin.x = (int)(xcoord - 3.0*res->pos.x);
	res->spin.y = (int)(ycoord - 3.0*res->pos.y);
	res->spin.z = (int)(zcoord - 3.0*res->pos.z);
      }else{
	cout<<"atom type not recognised: "<<atom<<endl;
	exit(1);
      }

    }else if(strcmp(type,"TER")==0){
      if(resnum>= LX ||resnum>= LY ||  resnum>= LZ){
	cout<<"warning: chain "<<cn<< " in file "<< s <<" may be too long for box "<<endl;
      }
      chains[cn]->N = resnum+1;
      
      chains[cn]->setFwdBkwdSpins();
     
      // chains[cn]->setRandomSpins();
      

      cn++;
      //cout <<"cn "<< cn<< endl;
      //cout<<line<<endl;
      chains[cn]=new Chain(this,cn);
    }else if(strcmp(type,"LIG")==0){
      chains[cn]->frozen = true;
      chains[cn]->locked = true;
    }else if(strcmp(type,"RIG")==0){
      chains[cn]->frozen = true;
    }
  }
  nChains = cn;
  applyPeriodicBoundaries();
  if(App::nativeStats){
    native = new Native(this);
  }
  stats.getLatticeStats(this);
  cout << "number of chains: "<< nChains<< endl;
}


void Lattice::resetStats(){
  stats.clean();
  stats.getLatticeStats(this);
}


void Lattice::applyPeriodicBoundaries(){
  for(int nc=0;nc<nChains;nc++){
    Chain * ch = chains[nc];
    for (int n=0;n<ch->N;n++){
      Residue * res =ch->residues[n];
      emptyLatticePos(res->pos);
      res->pos.periodicBoundary();  
      setResidue(res->pos,res);
    }
  }
}

void Lattice::writePDBMultiChain(string s){
  ofstream outs;
  cout<< "printing to "<< s<<endl;
  outs.open(s.c_str());
  writePDBMultiChain(outs,1);
  outs.close();
}

void Lattice::writePDBMultiChain(ofstream &outs,int moviestep){
  outs << "MODEL "<<moviestep<<endl;
  //outs <<  setw(7) << moviestep<<endl;
  for(int nc=0;nc<nChains;nc++){
    Chain * ch = chains[nc];
    if(ch->frozen && ch->locked){
      outs << "LIG"<<endl ;
    }else if(ch->frozen){
      outs << "RIG"<<endl;
    }
    for (int n=0;n<ch->N;n++){
      Residue * res =ch->residues[n];
      outs << "ATOM  " ;
      outs << setw(5)<< n;
      outs << "  CA  ";
      outs << setw(3)<< AA::int2aa[res->aa];
      outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
      outs << setw(4)<<n;
      outs << "   ";
      outs << setw(8)<< 3.0* (float) res->pos.x ;
      outs << setw(8)<< 3.0*(float) res->pos.y ;
      outs << setw(8)<< 3.0*(float) res->pos.z ;
      string occ = "1.00";
      string state = "11.00";
      if(res->state == inCoil){
	state = "22.00";
      }
      if(res->masked){
	occ="2.00";
      }
      outs<< "   "<<occ<<" "<<state;
      outs << endl;

      outs << "ATOM  " ;
      outs << setw(5)<< n;
      outs << "  CB  ";
      outs << setw(3)<< AA::int2aa[res->aa];
      outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
      outs << setw(4)<<n;
      outs << "   ";
      outs << setw(8)<< (3.0*(float) res->pos.x)+res->spin.x ;
      outs << setw(8)<< (3.0*(float) res->pos.y)+res->spin.y ;
      outs << setw(8)<< (3.0*(float) res->pos.z)+res->spin.z ;
      
      // take same state and mask as above
      outs<< "   "<<occ<<" "<<state;
      outs << endl;


    }
    outs <<"TER   " << endl; 
  }
  int printed=0;
  if(App::setAir){
    for (int x=0;x<LX;x++){
      for (int y=0;y<LY;y++){
	for (int z=0;z<LZ;z++){
	  if(r[x][y][z]==AIR && drand48()<0.03){
	    Residue * res=AIR;
	    outs << "ATOM  " ;
	    outs << setw(5)<< printed;
	    outs << "  CA  ";
	    outs << setw(3)<< "OOO";
	    outs << " "<< (char)'z';    // " A";
	    outs << setw(4)<<printed;
	    outs << "   ";
	    outs << setw(8)<< 3.0* (float)x ;
	    outs << setw(8)<< 3.0*(float) y ;
	    outs << setw(8)<< 3.0*(float) z ;
	    if(res->state == inCoil){
	      outs << "   1.00 22.00";
	    }else{
	      outs << "   1.00 11.00";
	    }
	    outs << endl;
	    printed++;
	  }
	}
      }
    }
  }
}

void Lattice::printPeriodicPDB(){

 ofstream outs;
 outs.open("periodic.pdb");
 outs << "MODEL "<<endl;
 for(int nc=0;nc<nChains;nc++){
   Chain * ch = chains[nc];
   int boundariesX[ch->N];
   bool leftX[ch->N];
   int boundariesY[ch->N];
   bool leftY[ch->N];
   int boundariesZ[ch->N];
   bool leftZ[ch->N];
   int nX=0;
   int nY=0;
   int nZ=0;
   //check for boundary crossings
   for (int n=0;n<ch->N;n++){
     Residue * res =ch->residues[n];
     cout<<n<<" "<<res->pos.x<<endl;
     if( res->pos.x==LX-1){
       if(n!=0 && ch->residues[n-1]->pos.x==0){
	 boundariesX[nX]=n; 
	 leftX[nX]=true;
	 nX++; 
	 cout<< "leftX at "<<n<<endl;
       }else if(n!= ch->N-1 && ch->residues[n+1]->pos.x==0){
	 boundariesX[nX]=n; 
	 leftX[nX]=false;
	 nX++; 
	 cout<< "rightX at "<<n<<endl;
       }
     }
     if( res->pos.y==LY-1){
       if(n!=0 && ch->residues[n-1]->pos.y==0){
	 boundariesY[nY]=n; 
	 leftY[nY]=true;
	 nY++; 
	 cout<< "leftY at "<<n<<endl;
       }else if(n!= ch->N-1 && ch->residues[n+1]->pos.y==0){
	 boundariesY[nY]=n; 
	 leftY[nY]=false;
	 nY++; 
	 cout<< "rightY at "<<n<<endl;
       }
     }
     if( res->pos.z==LZ-1){
       if(n!=0 && ch->residues[n-1]->pos.z==0){
	 boundariesZ[nZ]=n; 
	 leftZ[nZ]=true;
	 nZ++; 
	 cout<< "leftZ at "<<n<<endl;
       }else if(n!= ch->N-1 && ch->residues[n+1]->pos.z==0){
	 boundariesZ[nZ]=n; 
	 leftZ[nZ]=false;
	 nZ++; 
	 cout<< "rightZ at "<<n<<endl;
       }
     }
   }
   cout<<"checked boundaries"<<endl;
   int previousBoundaryX=-1;
   // add BOX length to left hanging regions
   for(int bc=0;bc<nX;bc++){
     if(leftX[bc]){
       cout<<"boundary at "<<boundariesX[bc] <<endl;
       for(int i = previousBoundaryX+1;i<boundariesX[bc];i++){
	 ch->residues[i]->pos.x += LX;
	 
       }
     }
     previousBoundaryX = boundariesX[bc];
   }
   if(nX>0 &&  !leftX[nX-1]){
     cout<<"add at end part "<<endl;
     for(int i = previousBoundaryX+1;i<ch->N;i++){
       ch->residues[i]->pos.x += LX;
     }
   }
   int previousBoundaryY=-1;
   // add BOX length to left hanging regions
   for(int bc=0;bc<nY;bc++){
     if(leftY[bc]){
       cout<<"boundary at "<<boundariesY[bc] <<endl;
       for(int i = previousBoundaryY+1;i<boundariesY[bc];i++){
	 ch->residues[i]->pos.y += LY;
	 
       }
     }
     previousBoundaryY = boundariesY[bc];
   }
   if(nY>0 &&  !leftY[nY-1]){
     cout<<"add at end part "<<endl;
     for(int i = previousBoundaryY+1;i<ch->N;i++){
       ch->residues[i]->pos.y += LY;
     }
   }
   int previousBoundaryZ=-1;
   // add BOX length to left hanging regions
   for(int bc=0;bc<nZ;bc++){
     if(leftZ[bc]){
       cout<<"boundary at "<<boundariesZ[bc] <<endl;
       for(int i = previousBoundaryZ+1;i<boundariesZ[bc];i++){
	 ch->residues[i]->pos.z += LZ;
	 
       }
     }
     previousBoundaryZ = boundariesZ[bc];
   }
   if(nZ>0 &&  !leftZ[nZ-1]){
     cout<<"add at end part "<<endl;
     for(int i = previousBoundaryZ+1;i<ch->N;i++){
       ch->residues[i]->pos.z += LZ;
     }
   }
   // print chain
   for (int n=0;n<ch->N;n++){
     Residue * res =ch->residues[n];
     outs << "ATOM  " ;
     outs << setw(5)<< n;
     outs << "  CA  ";
     outs << setw(3)<< AA::int2aa[res->aa];
     outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
     outs << setw(4)<<n;
     outs << "   ";
     outs << setw(8)<< (float) res->pos.x ;
     outs << setw(8)<< (float) res->pos.y ;
     outs << setw(8)<< (float) res->pos.z ;
     outs << "   1.00 22.00";
     outs << endl;
   }
   outs <<"TER   " << endl;
    for (int n=0;n<ch->N;n++){
     Residue * res =ch->residues[n];
     outs << "ATOM  " ;
     outs << setw(5)<< n;
     outs << "  CA  ";
     outs << setw(3)<< AA::int2aa[res->aa];
     outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
     outs << setw(4)<<n;
      outs << "   ";
      outs << setw(8)<< (float) res->pos.x +LX;
      outs << setw(8)<< (float) res->pos.y ;
      outs << setw(8)<< (float) res->pos.z ;
      outs << "   1.00 22.00";
      outs << endl;
   }
outs <<"TER   " << endl; 
for (int n=0;n<ch->N;n++){
     Residue * res =ch->residues[n];
     outs << "ATOM  " ;
     outs << setw(5)<< n;
     outs << "  CA  ";
     outs << setw(3)<< AA::int2aa[res->aa];
      outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
      outs << setw(4)<<n;
      outs << "   ";
      outs << setw(8)<< (float) res->pos.x ;
      outs << setw(8)<< (float) res->pos.y +LY;
      outs << setw(8)<< (float) res->pos.z ;
      outs << "   1.00 22.00";
      outs << endl;
   }
outs <<"TER   " << endl; 

for (int n=0;n<ch->N;n++){
     Residue * res =ch->residues[n];
     outs << "ATOM  " ;
     outs << setw(5)<< n;
     outs << "  CA  ";
     outs << setw(3)<< AA::int2aa[res->aa];
      outs << " "<< (char)('A' +(res->chainNum)%('A' - 'z'));    // " A";
      outs << setw(4)<<n;
      outs << "   ";
      outs << setw(8)<< (float) res->pos.x ;
      outs << setw(8)<< (float) res->pos.y ;
      outs << setw(8)<< (float) res->pos.z+LZ ;
      outs << "   1.00 22.00";
      outs << endl;
      }

   outs <<"TER   " << endl; 
 }
 outs.close();
}




int Lattice::insertChain(Config * config, int * sequence,double beta){
  //cout<<"inserting chain at: "<< nChains <<endl;
  int chainNum = nChains;
  if(chainNum == MAX_CHAINS){cout<< "ERROR too many chains"<<endl;exit(0);}
  Chain * chain = new Chain (this,chainNum);  
  //chain->setBeta(beta);
  //cout <<"BETA "<<beta<<endl; 
  chain->frozen = false;
  chain->locked = false;
  chain->N = config->nTotal;
  for(int n=0;n<config->nTotal;n++){
    Residue * res = new Residue;
    res->pos.x = config->positions[n].x;
    res->pos.y = config->positions[n].y;
    res->pos.z = config->positions[n].z;
    res->spin=config->spins[n];
    res->state=config->states[n];
    res->aa = sequence[n];
    res->chainNum = chainNum;
    res->n =n;
      
    chain->residues[n]=res;
    r[res->pos.x][res->pos.y][res->pos.z] = chain->residues[n];
  }
  chains[chainNum]=chain;
  nChains++;
  chain->setFwdBkwdSpins();
  Stats s;
  s.get_Eint_Cint(chains[chainNum],this);
  stats += s;
  //freeChains++;
  return 0;
}

int Lattice::deleteChain(int chain){
  //cout<<"deleting chain: "<< chain<<endl;
  // delete residues from lattice

  Stats s;
  s.get_Eint_Cint(chains[chain],this);
  stats -= s;

  for(int n=0;n<chains[chain]->N;n++){
    Residue * res=chains[chain]->residues[n];
    r[res->pos.x][res->pos.y][res->pos.z]=NULL;
  }
  delete chains[chain];
  nChains--;
  //freeChains--;

  // replace last chain, to new position, changeing chainNum in residues
  if(chain != nChains){
    for(int n=0;n<chains[nChains]->N;n++){
      chains[nChains]->residues[n]->chainNum=chain;
    }
    chains[chain] =chains[nChains];
    chains[chain]->setChainNum(chain);
  }
  //cout<<"finished deleting chain: "<< chain<<endl;
  return 0;
}






bool Lattice::checkStats(){
  Stats newStats;
  newStats.getLatticeStats(this);
  if(stats != newStats){
    cout<< "ERROR IN STATS"<<endl;
    cout<< "calculated stats:"<<endl;
    newStats.printCout();
    cout<< "additive stats:"<<endl;
    stats.printCout();
    writePDBMultiChain("errorInStats.pdb");
    exit(1);
    return false;
  }else{
    return true;
  }

}


void Lattice::oldReadPDBMultiChain(string s){
  ifstream infile(s.c_str());
  if (!infile.is_open()){
    cout << "can't open file "<< s <<endl;
    exit(0);
  }
  int resnum;
  int cn =0;
  char line[500];
  char type[5],restype[4], atom[4];;
  float xcoord,ycoord,zcoord;
  float occupancy;
  // bool lig=false;
  chains[cn]=new Chain(this,cn);
  
  while(!infile.eof()){
    if(cn>= MAX_CHAINS){
      cout<<"Too many chains in "<< s <<endl;
      exit(0);
    }
    infile.getline(line,2000);
    type[0]=' ';
    sscanf(line,"%4s%*9c%2s%*2c%3s%*2c%4d%*4c%9f%9f%9f%4f\n",type,atom,restype,&resnum,&xcoord,&ycoord,&zcoord,&occupancy);
    if ((strcmp(type,"ATOM")==0 && strcmp(atom,"CA")==0) &&(strcmp(line,"")!=0)){

      if(resnum>= MAX_RES){
	cout<<"Too many residues in chain "<<cn<< " in file "<< s <<endl;
	exit(0);
      }
      
      chains[cn]->residues[resnum] = new Residue;
      chains[cn]->residues[resnum]->pos.x =(int) xcoord;
      chains[cn]->residues[resnum]->pos.y = (int) ycoord;
      chains[cn]->residues[resnum]->pos.z =(int) zcoord;
      chains[cn]->residues[resnum]->aa = aaInt->stringtoAA(restype);
      chains[cn]->residues[resnum]->n =resnum;
      chains[cn]->residues[resnum]->chainNum = cn;
      r[(int)xcoord][(int)ycoord][(int)zcoord]= chains[cn]->residues[resnum];
    }else if(strcmp(type,"TER")==0){
      if(resnum>= LX ||resnum>= LY ||  resnum>= LZ){
	cout<<"warning: chain "<<cn<< " in file "<< s <<" may be too long for box "<<endl;
      }
      chains[cn]->N = resnum+1;
      
      chains[cn]->setFwdBkwdSpins();
     
      chains[cn]->setRandomSpins();
      

      cn++;
      //cout <<"cn "<< cn<< endl;
      //cout<<line<<endl;
      chains[cn]=new Chain(this,cn);
    }else if(strcmp(type,"LIG")==0){
      chains[cn]->frozen = true;
      chains[cn]->locked = true;
    }else if(strcmp(type,"RIG")==0){
      chains[cn]->frozen = true;
    }
  }
  nChains = cn;

  stats.getLatticeStats(this);
  cout << "number of chains: "<< nChains<< endl;
}



void Lattice::setNative(string fn_native){
  if(App::nativeStats){
  cout<< "SETTING NATIVE"<<endl;
  Lattice * tmp = new Lattice;
  tmp->setAA(App::fn_aa, App::indx_water);
  tmp->readPDBMultiChain(fn_native);
  native->setNative(tmp);
  delete tmp;  // TO DO delete properly
  stats.getLatticeStats(this);
  }
}

