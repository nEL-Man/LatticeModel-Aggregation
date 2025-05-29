#include "Pos.hpp"


#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>


using namespace std;


void Pos::printCout(){
  cout<<"x:"<<x<<" y:"<<y<<" z:"<<z<<endl;
}


string Pos::toString(){
  std::stringstream myStream;
  myStream<<"x:"<<x<<" y:"<<y<<" z:"<<z;
  return myStream.str();
}
