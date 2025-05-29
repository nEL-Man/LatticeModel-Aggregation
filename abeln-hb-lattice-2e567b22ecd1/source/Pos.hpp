/**
 * \file   Pos.hpp
 * \date   May 2013
 * \author Sanne Abeln
 * \brief  Defines class Pos
 */



#ifndef _POS_HPP
#define _POS_HPP

#include <stdio.h>

#define XDIR 0 
#define YDIR 1 
#define ZDIR 2 

#define LX 30
#define LY 30
#define LZ 30

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string.h>


using namespace std;


/// Class storing position vectors
/// its methods take care of the periodic boundaries
class Pos{
public:
  /// Constructor, taking coordinates
  Pos(int tx, int ty, int tz );
  /// deault constructor
  Pos();
  Pos operator+	(const Pos & v) const;
  Pos operator-	(const Pos & v) const;
  Pos operator- ()const;
  Pos& operator=  (const Pos & p);
  bool operator!= (const Pos &p)const;
  bool operator== (const Pos &p)const;
  void periodicBoundary();
  void periodicSubtraction(Pos p1,Pos p2);
  void printCout();
  string toString();
  static bool orthogonal(Pos p1,Pos p2);
  const int operator[] ( const int indx )const;
  int& operator[] ( const int indx );
  // --- data structure ---
  union
  {
    struct
    {
      int x;
      int y;
      int z;
    };
    int	xyz[3];
  };
private:
  void copy(const Pos& p);
};


// Non-member functions
int getAngleUnitV(Pos p1,Pos p2);
int getAnglePositionV(Pos p1,Pos p2,Pos p3);
Pos pickNeighbour90(Pos p1,Pos p3, int randomN);
Pos pickNeighbour120(Pos p1,Pos p2,Pos p3);

inline
Pos::Pos(int tx,int ty,int tz){
  x=tx;y=ty;z=tz;
}

inline
Pos::Pos(){
  x=0;y=0;z=0;
}



inline
Pos Pos::operator+ (const Pos & v)const{
  return Pos(x+v.x,y+v.y,z+v.z);
}

inline
Pos Pos::operator- (const Pos & v)const{
  return Pos(x-v.x,y-v.y,z-v.z);
}

inline
Pos Pos::operator-()const{
  return Pos(-x,-y,-z);
}


inline
bool Pos::operator!= (const Pos &p)const{
  return (x!=p.x || y!=p.y || z!=p.z);
}

inline
bool Pos::operator== (const Pos &p)const{
  return (x==p.x && y==p.y && z==p.z);
}



inline
const int Pos::operator[] (const int indx)const{
  return xyz[indx];
}

inline
int &  Pos::operator[] (const int indx){
  return xyz[indx];
}



inline
Pos&	Pos::operator= 	( const Pos & p )
{
	copy( p );
	return (*this);
}

inline
void Pos::copy ( const Pos & p ){
	x = p.x;
	y = p.y;
	z = p.z;
}





inline
void Pos::periodicBoundary(){
  if(x>=LX){
    x =x-LX;
  }
  if(x<0){
    x =LX + x;
  }
 if(y>=LY){
    y =y-LY;
  }
  if(y<0){
    y =LY + y;
  }
  if(z>=LZ){
    z =z-LZ;
  }
  if(z<0){
    z =LZ + z;
  }
}

inline
void Pos::periodicSubtraction(Pos p1,Pos p2){
  x=p1.x-p2.x;
  y=p1.y-p2.y;
  z=p1.z-p2.z;
  if(x > (LX/2)){x = x -LX;}
  else if(x< -(LX/2)){x=x+LX;}
  if(y > (LY/2)){y=y -LY;}
  else if(y< -(LY/2)){y=y+LY;}
  if(z > (LZ/2)){z=z -LZ;}
  else if(z< -(LZ/2)){z=z+LZ;}
}


inline 
bool Pos::orthogonal(Pos p1, Pos p2){
  bool dot_product = p1.x*p2.x && p1.y*p2.y && p1.z*p2.z;
  return (dot_product ==0);
}



#endif
