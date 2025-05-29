/**
 * \file   ptemp.cpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines methods for parallel tempering using openmpi also contains the "main" function
 */


#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cstdlib>


#include "GranCan.hpp"
#include "App.hpp"
#include "ptemp.h"
#include "Design.hpp"
#include "MonteCarlo.hpp"


#define LISA

using namespace std;
//using namespace __gnu_cxx;

//mpiCC tstmpi.cpp -o tstmpi



// MESSAGE TYPES
#define SWAP_REQUEST 1
#define ALL_DONE 2



// FORWARD DECLARATIONS

bool trySwap(int,int);
int initBetas(double minBeta,double maxBeta);

// GLOBALS
int MASTER =0;
int myRank;
bool LAMBDA_COR = false;


double betas[MAX_PROCS];
int betaID2node[MAX_PROCS];
int node2betaID[MAX_PROCS];
int energies[MAX_PROCS];
int nChains[MAX_PROCS];
int swapacc[MAX_PROCS]={0};
int swaptry[MAX_PROCS]={0};


double myBeta;
int num_procs;

bool swapT=true;

//Rates * transitionRates;

//double fraction_swaps=0.01;


int main(int argc, char *argv[]) {

  ///////////////////
  // contains instances of classes:
  // Design
  Design * design=NULL;
  
  // MonteCarlo
  MonteCarlo * mcSim = NULL;

  // App is used to store global variables and command line arguments
  //
  ////////////////////

  //int steps=0;
  char name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI::Get_processor_name(name, namelen);
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /// get process rank
  myRank = MPI::COMM_WORLD.Get_rank ( );
  App::myRank= myRank;

  /// check number of processes
  if(num_procs > MAX_PROCS){
    cout<< "too many processes: "<<num_procs;
    exit(1);
  }
  printf("Process %d out of %d realp %s \n", myRank, num_procs,name);

  /// initialise random number generator with rank
  srand48(myRank);  

  /// pass arguments to App
  App::init(argc, argv);


  /// initialise Monte Carlo and/or  Design
  if(App::designProgram){
    design= new Design();
    if(! App::useOldFormat){
      design->init(App::fn_in);
    }else{
      design->initOldFormat(App::fn_in);
    }
  }
 
  
  mcSim = App::getMonteCarlo();

  cout<< "finished intilialising"<< endl;
  if(myRank==MASTER){
    initBetas(App::minBeta,App::maxBeta);
  }

  MPI::COMM_WORLD.Bcast(&node2betaID,num_procs,MPI::INT, MASTER);
  MPI::COMM_WORLD.Bcast(&betas,num_procs,MPI::DOUBLE, MASTER);
  MPI::COMM_WORLD.Bcast(&betaID2node,num_procs,MPI::INT, MASTER);
  // MPI::COMM_WORLD.Barrier();
  myBeta=betas[node2betaID[myRank]];

  cout<< "betas set"<< endl;

  for(int swaps=0;swaps<mcSim->total_swaps;swaps++){
    //cout<< "MC process: "<<myRank<< " at beta: "<< myBeta<<" ";
    //cout<< "betaID "<< node2betaID[myRank]<<endl;  
    
    //should check if beta is changed ... 
    bool myBetaChanged = true;
    
    int myEnergy;
    if(App::designProgram){    
      design->designProcedure(myBeta,mcSim->steps);
    }else{
      myEnergy = mcSim->MC( myBeta,node2betaID[myRank],myBetaChanged);
    }


    MPI::COMM_WORLD.Gather(&myEnergy,1,MPI::INT,
			   &energies,1, MPI::INT,
			   MASTER); 
    
    int myNChains = mcSim->grancan->getNumFreeChains();
    MPI::COMM_WORLD.Gather(&myNChains,1,MPI::INT,
			   &nChains,1, MPI::INT,
			   MASTER);
    if(swapT){
      if(myRank==MASTER){
	//cout <<endl;
	int start =0;// swap even
	if (drand48()<0.5) start=1; //swap uneven
	for(int i=start;i<num_procs-1;i=i+2){
	  trySwap(i,i+1);
	}
      }
    }
    MPI::COMM_WORLD.Bcast(&node2betaID,num_procs,MPI::INT, MASTER);
    MPI::COMM_WORLD.Bcast(&betaID2node,num_procs,MPI::INT, MASTER);
    myBeta=betas[node2betaID[myRank]];
  }
  //double myEnergy = grancan->MC(steps,myBeta,node2betaID[myRank],true);

  /// finalize
  if(App::designProgram){    
    design->finalize();
  }else{
    mcSim->getMCStats(num_procs,myRank);
  }
  if(myRank==MASTER){
    for(int  i = 0 ; i<num_procs;i++) {
      cout << "beta " << i << "trials=" << swaptry[i] <<" acc=" << swapacc[i] << " ratio=" << (double)swapacc[i]/swaptry[i]  << endl;
    }
  }

 
#ifdef LISA
  //  cout << "waiting for barrier rank "<<myRank<<endl;
  //MPI::COMM_WORLD.Barrier();
  cout << "finalizing rank "<<myRank<<endl;
  MPI::Finalize();
#endif
  cout << "close down process "<<myRank<<endl;
  exit(0);
}


//double zz = log(eBetaMu);

bool trySwap(int betaID1, int betaID2){
  cout <<"trial: "<<betaID1<<" "<<betaID2<<endl;
  int node1 = betaID2node[betaID1];
  int node2 = betaID2node[betaID2];
  double beta1 = betas[betaID1];
  double beta2 = betas[betaID2];
  double realBeta1=100*beta1;
  double realBeta2=100*beta2;

  double energy1 =energies[node1];
  double energy2 =energies[node2];

  int nChains1 = nChains[node1];
  int nChains2 = nChains[node2];


  // cout<< "energies: " <<energy1<<" "<<energy2<<endl;
  double dE = (double) (energy1-energy2);
  double dB = beta1-beta2;
  //double pAccSwap = exp(zz*dN + dE*dB);
  double prefactor;  
  if(LAMBDA_COR){
    prefactor = pow(realBeta1/realBeta2,(3/2)*(nChains1-nChains2));
  }else{
    prefactor=1.0;
  }

  double pAccSwap = prefactor*exp( dE*dB);
  swaptry[betaID1]++;
  swaptry[betaID2]++;

  // cout << "TRY SWAP between " <<node1 << " and " << node2<<endl;
  if(drand48()<pAccSwap){
    cout << "SWAPPING " << beta1 <<" with energy: "<< energy1 <<" at node "<<node1;
    cout << " with " << beta2<<" with energy: "<< energy2<< " at node " <<node2<<endl;
    betaID2node[betaID1] =node2;
    betaID2node[betaID2] =node1;
    node2betaID[node1]=betaID2;
    node2betaID[node2]=betaID1;
    swapacc[betaID1]++;
    swapacc[betaID2]++;

    //transitionRates->changeBeta(Cnat,ligC);
    return true;
  }else{
    return false;
  }
}




int initBetas(double minBeta,double maxBeta){
  if(num_procs ==1){
    betas[0]= (minBeta + maxBeta)/2;
     betaID2node[0]=0;
     node2betaID[0]=0;
  }else{
    if(!App::linear_T){
      double step=(double)(maxBeta - minBeta)/(double) (num_procs-1);
      for (int i=0;i<num_procs;i++){
	betas[i]= minBeta + i*step;
	betaID2node[i]=i;
	node2betaID[i]=i;
	cout<< "node "<<i<<" beta "<<betas[i]<<" T "<<0.01/ betas[i]<<endl;
      }
    }else{
      double maxT =0.01/minBeta;
      double minT =0.01/maxBeta;
      double temperatures[MAX_PROCS];
      double step=(double)(maxT - minT)/(double) (num_procs-1);
      for (int i=0;i<num_procs;i++){
	temperatures[num_procs - i-1]= minT + i*step;
      }
      for (int i=0;i<num_procs;i++){
	betas[i]= 0.01/temperatures[i];
	betaID2node[i]=i;
	node2betaID[i]=i;
	cout<< "node "<<i<<" beta "<<betas[i]<<" T "<<temperatures[i]<<endl;
      }
    }
  }
  return 0;
}
