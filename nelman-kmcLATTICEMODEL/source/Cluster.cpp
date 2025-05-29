#include "Cluster.hpp"
#include "Lattice.hpp"
#include "Residue.hpp"
#include "Chain.hpp"
#include "AA.hpp"
#include "Stats.hpp"
#include "StatsMethods.hpp"
#include "App.hpp"


#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


Cluster::Cluster(){
}



int Cluster::oldClusterSize=0;

void Cluster::reset(){
  largestCluster=0;
  totalClusters=0;
  totalVisited=0;
  for(int i=0 ;i<MAX_CHAINS;i++){
    visited[i] = -1;
    clustNum[i]= -1;
  }
  for(int i=0 ;i<MAX_CHAINS;i++){
    clusterCount[i] = 0;
  }
}


void Cluster::startCounting(Lattice * l1){
  l=l1;
  reset();
  int curCluster=0;
  // cout<<"starting cluster count!!"<<endl;
  totalChains = l->nChains;
  for(int i=l->nChains -1;i>=0;i--){
    if(visited[i]<0){
      visit(i,curCluster);
      if(clusterCount[curCluster]>clusterCount[largestCluster]){
	largestCluster=curCluster;
      }
      curCluster++;
    }
  }
  totalClusters = curCluster;
}


void  Cluster::visit(int chain, int cluster){
  clusterCount[cluster]++;
  totalVisited++;
  visited[chain]=cluster;
  clustNum[chain]=cluster;
  for(int n=0; n< l->chains[chain]->N;n++){
    if(totalVisited >= l->nChains){
      //cout<<"early return"<<endl;
      return;
    }
    Residue * res= l->chains[chain]->residues[n];
    //if(res->aa !=13){
    for(int k=0;k<6;k++){
      Pos nb =local[k]+ res->pos;
      nb.periodicBoundary();
      Residue * neighbour =l->getResidue(nb); 
      //if(neighbour!=NULL && neighbour->aa != 13){
      if(neighbour!=NULL && neighbour != l->AIR){
	if(!res->masked &&  !neighbour->masked){
	  int nChain = neighbour->chainNum;
	  if(visited[nChain]<0){
	    visit(nChain,cluster);
	  }
	}
      }
    }
#ifdef CROSSNBS
    Pos * crossNbs = Stats::getCrossNbsSpinDir(res->spin);
    for(int k=0;k<4;k++){
      Pos posNB =crossNbs[k]+ res->pos;
      posNB.periodicBoundary();
      Residue * resNB = l->getResidue(posNB);
      if(resNB!=NULL  && neighbour != l->AIR){
	bool sameChain = resNB->chainNum==res->chainNum;
	if(!sameChain){
	  bool spinsOpposite= (resNB->spin == (Pos(0,0,0)- res->spin));
	  if(spinsOpposite){
	    int nChain = resNB->chainNum;
	    if(visited[nChain]<0){
	      visit(nChain,cluster);
	    }
	  }
	}
      }
    }
    //delete crossNbs;
#endif
  }
}



void Cluster::printCout(){
  for(int c=0;c<totalClusters;c++){
    cout <<"cluster "<<c<< " contains ";
    cout << clusterCount[c] << " chains"<<endl;
  }
}


ClusterInfo Cluster::checkCluster(Lattice *l1,double beta){
  l=l1;
  startCounting(l);
  ClusterInfo ci;
  int lowLim = (int)(0.95* (double)oldClusterSize);
  int highLim = (int)(1.05* (double)oldClusterSize);

  if(clusterCount[largestCluster] < lowLim || clusterCount[largestCluster]> highLim ){
    cout <<"Largest cluster "<< largestCluster << " size "<<clusterCount[largestCluster]<<endl; 
    if(App::writeCluster){
      writePDBClust(largestCluster,beta);
    }
    oldClusterSize= clusterCount[largestCluster];
  }
  //printCout();
  ci.clusters=totalClusters;
  ci.largestCluster = clusterCount[largestCluster];
  ci.totalChains = l->nChains;
  return ci;
}



void Cluster::writePDBClust(int clust, double beta){
  std::stringstream buf;
  buf << std::fixed << std::setprecision(2) << 0.01/beta; 
  string stT= buf.str();
  int size=clusterCount[clust]; 

  std::stringstream outPDB;
  outPDB << "clusT"<<stT<<"S"<<size<<"all.pdb";
  l->writePDBMultiChain(outPDB.str());

  std::stringstream outPDBclust;
  outPDBclust << "clusT"<<stT<<"S"<<size<<".pdb";
  ofstream outs;
  cout<< "printing to "<< outPDBclust.str()<<endl;
  outs.open(outPDBclust.str().c_str()); 
  outs << "MODEL "<<endl;
  //outs <<  setw(7) << moviestep<<endl;
  for(int nc=0;nc<l->nChains;nc++){
    if(visited[nc]==clust){
      Chain * ch = l->chains[nc];
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
	if(res->state == inCoil){
	  outs << "   1.00 22.00";
	}else{
	  outs << "   1.00 11.00";
	}
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
	if(res->state == inCoil){
	  outs << "   1.00 22.00";
	}else{
	  outs << "   1.00 11.00";
	}
	outs << endl;
      }
      outs <<"TER   " << endl; 
    }
  }
  outs.close();

  
  std::stringstream outClust;
  outClust << "clusT"<<stT<<"S"<<size<<".out";
  ofstream outfile(outClust.str().c_str());
  outfile << "beta: "<<beta<<endl;
  outfile << "beta: "<<stT<<endl;
  outfile << "eBetaMu: "<<App::eBetaMu<<endl;
  outfile << "LX LY LZ: "<<LX<<" "<<LY<<" "<<LZ<<endl;
  //print stats
  outfile << l->stats.print2string()<<endl;
  // print all stats
  for(int cc=0;cc<totalClusters;cc++){
    outfile <<"cluster "<<cc<< " contains ";
    outfile << clusterCount[cc] << " chains"<<endl;
  }
  outfile.close();
}


vector<int> Cluster::getChains(int cluster_number){
  vector<int> chainList; 
  for(int nc=0;nc<totalChains;nc++){
     if(cluster_number == clustNum[nc]){
       chainList.push_back(nc);
     }
   }
  return chainList;
}
