#include "AA.hpp"

#include <math.h>
#include <fstream>
#include <vector>
#include <sstream>

vector<string> split(const string& str, const string& delimiters = " ");
string toUpper(string in);

int averageInt[MAXAA];

int AA::WATER=0;
int AA::NUMAA=0;
int AA::designWater=0;

string AA::int2aa[MAXAA];



AA::AA(string filename, int indx_water){
  
  WATER = indx_water; 
 
  ifstream aaFile (filename.c_str());
  if (!aaFile.is_open()){
    cout << "could not open file: "<<filename.c_str()<<  endl;
    exit(1);
  }
  int i=0;
  while(!aaFile.eof()){
    string line;
    getline(aaFile,line);
    //split line into fields
    vector<string> fields = split(line, " ");
    //check if line not empty
    if(fields.size() > 0){
      //get amino acid
      string aa = fields[0];
      aa=toUpper(aa);
      //store amino acid
      int2aa[i]= aa;
      if(aa=="HOH") designWater=i;
      if(fields.size() > MAXAA +1){
	cout<<"too many amino acids in file (1) "<<filename<<" "<<fields.size() <<endl;
	cout<<line<<endl;
	cout<<fields[fields.size()-1]<<endl;
	exit(1);
      }
      //loop through all values matrix values
      for(int j=i+1;j< (int) fields.size();j++){
	//convert string to int, value * 100   
	std::istringstream is(fields[j]);
	double val ; is>>val;
	aaInt[i][j-1]=(int) rint(100 * val);
	aaInt[j-1][i]=(int) rint(100 * val);
      }
      i++;
    }
  }
  NUMAA=i;
  cout<<NUMAA<< " amino acids"<<endl;
  if(NUMAA > MAXAA){
    cout<<"too many amino acids in file (2) "<<filename<<" "<<NUMAA <<endl;
    exit(1);
  }

  if(indx_water >= NUMAA){
    cout<<"-indxWater "<<indx_water<<" exceeds range in: ";
    cout<< filename<<" with "<<NUMAA<<" amino acids "<<endl;
    exit(1);
  }
  
#ifdef DEBUG
  for(int i=0;i<NUMAA;i++){
    cout << int2aa[i] <<"\t";
    for(int j=0;j<NUMAA;j++){
      cout<< aaInt[i][j]<<"\t";
    }
    cout<<endl;
  }
#endif
 // OPT would need reordering with other matrix input
  aaFile.close();
}


int AA::stringtoAA(string s){
  s=toUpper(s);
  for(int i=0;i<NUMAA;i++){
    if(!s.compare(int2aa[i])){
      return i;
    }
  }
  cout << "aa not found:" <<s<<endl;
  exit(1);
  return(0);
}

vector<string> split(const string& str, const string& delimiters){
  vector<string> tokens;
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  return tokens;
}


string toUpper(string in){
  string out=in;
  for(int i=0;i<(int) in.size();i++){
    out[i]=toupper(in[i]);
  }
  return out;
}
