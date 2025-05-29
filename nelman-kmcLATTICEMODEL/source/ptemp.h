/**
 * \file   ptemp.h
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines methods for parallel tempering using openmpi also contains the "main" function
 */


#ifndef _ptemp_h_
#define _ptemp_h_


#define MAX_PROCS 32


extern double betas[MAX_PROCS];
extern int myRank;
extern bool swapT;
extern bool LAMBDA_COR;

#endif
