/**
 * \file   Cluster.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Cluster
 */

#ifndef _Cluster_H_
#define _Cluster_H_

#include <vector>

#include "Lattice.hpp"
#include "EnergyMap.hpp"


using namespace std;

/// Class containing information over the size of all clusters formed in a GranCan simulation 
class Cluster{
public:
  Cluster();
  void startCounting(Lattice * l1);
  void printCout();
  ClusterInfo checkCluster(Lattice *l,double beta);
  int getTotalClusters();

  vector<int> getChains(int cluster_number);
  // int checkCluster(Lattice *l,double beta, int required_cluster_size);
private:
  void writePDBClust(int clust, double beta);
  void visit(int chain, int cluster);
  void reset();
  int totalChains;
  int clusterCount[MAX_CHAINS];
  int totalClusters;
  int visited[MAX_CHAINS];
  int clustNum[MAX_CHAINS];
  int totalVisited;
  static  int oldClusterSize;
  int largestCluster;
  Lattice *l;
};

inline
int Cluster::getTotalClusters(){
  return totalClusters;
}



#endif
