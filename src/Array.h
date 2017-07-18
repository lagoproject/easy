//----------------------------------------------------------------------
/*
  File Array.h 
 */
//----------------------------------------------------------------------

#ifndef ARRAY_H
#define ARRAY_H

#include <cmath>
#include <vector>
#include <string>
#include "TObject.h" 
#include "Constants.h" 
#include "Station.h" 

//----------------------------------------------------------------
/*
  class Array
  The array of Stations used for the simulation.
  An array can be the full array or a preselected array of station.
  An array can be created 3 different ways :
      - by reading a list of stations in a file 
      - a default array can be created automatically ( then the 
        argument is the number of stations wanted in the array )
      - from a vector of stations 
*/
//----------------------------------------------------------------
    
class Array : public TObject
{
 private:
  
 public:
  
  vector<Station> fStationList;
  Double_t fNorMin;
  Double_t fNorMax;
  Double_t fEasMin;
  Double_t fEasMax;
  Double_t fNorMean;
  Double_t fEasMean;
  Double_t fNorWidth;
  Double_t fEasWidth;
 
  Array(){};
  Array(vector<Station> tanklist);
  Array(Int_t size);
  Array(Double_t theta,Double_t phi,Double_t E0, Double_t N0);
  Array(string arrayfilename);
  ~Array(){};
  void Fill(Station* sta);
  Int_t Size();
 
  ClassDef(Array,1)
    
};
extern Array *gArray;
extern Array *gHexagArray;
#endif
