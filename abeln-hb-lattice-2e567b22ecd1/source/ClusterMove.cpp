#define MAX_TRANSLATION 10


int ClusterMove(){
  /// calculate clusters
  Cluster cl;
  cl.startCounting(mainLattice);
  int numClusts = cl.getTotalClusters();

  /// pick random clusters
  int cluster2move = (int) floor(numClusts * drand48());

  /// get all chains in cluster
  vector<int> chainList = cl.getChains(cluster2move);

  /// set rotation and translation
  //int dir_translation= (int)floor(6 * drand48());
  //int magnitude_translation = (int)floor(MAX_TRANSLATION * drand48());
  int rotationalDir= (int)floor(6 * drand48());

  /// choose random pivot on lattice
  int xpivot =  (int)floor(LX * drand48()); 
  int ypivot =  (int)floor(LY * drand48());
  int zpivot =  (int)floor(LZ * drand48());
  Pos posPivot =  Pos(xpivot,ypivot,zpivot);

  vector< vector Pos>> newPos;
  
  /// for each chain, perform rotation and translation,
  /// while keeping old coordnates
  /// check for clashes or touches: reject, restore old coordinates and sample
  for(int i=0; i< chainList.size();i++){
    int cn = chainList[i];
    Chain * ch = l->chains[cn];    
    for(int n=ch->startRes; n<=ch->endRes;n++){
      Residue * res = ch->residues[n];
      posNew[cn][n]=  ch->rotationPosition(res->pos,rotationalDir,posPivot);
      if(l->checkForClashesAndTouches(posNew[cn][n])){
	/// reject move
	l->sampleCurrentStats();
	return -1;
      }      

    }
  }
  /// now we know there are no clashes our touches
  /// and the move can be performed
  /// note that the move will alsways get accepted, since dE=0


  /// Perform move, update lattice
  for(int i=0; i< chainList.size();i++){
    int cn = chainList[i];
    Chain * ch = l->chains[cn];
    for(int n=ch->startRes; n<=ch->endRes;n++){
      Residue * res = ch->residues[n];
      Pos newSpin = ch->rotationSpin(res->spin,rotationalDir);
      l->emptyLatticePos(res->pos);

      res->pos = posNew[cn][n];
      res->spin = newSpin;
      l->setResidue(posNew[cn][n],res)
    }
  }
}


