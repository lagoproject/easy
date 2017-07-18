/*
  Array.cc
  implementation file for class Array
*/

#include "Array.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <strstream>
#include <stdlib.h>
#include "EasySimConfig.h"
//#include "Utils.h"

ClassImp(Array)

Array *gArray;
Array *gHexagArray;

Array::Array(vector<Station> tanklist)
{
  Double_t north,east;
  fStationList= tanklist;
  fNorMin=99999999;
  fNorMax=0;
  fEasMin=99999999;
  fEasMax=0;
  Int_t size = fStationList.size();
  for (Int_t i=0;i<size;i++){
    north = fStationList[i].fNorthing;
    east = fStationList[i].fEasting;
    if(north>fNorMax)fNorMax = north;
    if(north<fNorMin)fNorMin = north;  
    if(east>fEasMax)fEasMax = east;
    if(east<fEasMin)fEasMin = east;     
  }
  fNorMean = (fNorMax+fNorMin)/2. ;
  fEasMean = (fEasMax+fEasMin)/2. ;
  fNorWidth = fNorMax-fNorMin ;
  fEasWidth = fEasMax-fEasMin ;
}




Array::Array(Int_t size)
{
  fNorMin=99999999;
  fNorMax=0;
  fEasMin=99999999;
  fEasMax=0;
   Int_t k = 101;
  Int_t NMIN = 6077556;
  Int_t EMIN = 455129;
  Int_t n;
   ofstream stafile;
  stafile.open("stationsdefault.txt",ios::out);
  if(theConfig ()->ArrayMode == "DEFAULTSOUTH"){
    
    n=(Int_t)sqrt(2.*size);
    
    // cout<<"Write a default array of "<<size<<" stations"<<endl;
    for(Int_t i=0;i<n;i++)
      {
      for(Int_t j=0;j<n;j+=2) 
	{
	  Station * sta= new Station (k,NMIN+i*1300,EMIN+j*750+750*(i%2),1400);
	  stafile<<k<<"\t"<<NMIN+i*1300<<"\t"<<EMIN+j*750+750*(i%2)<<"\t 1400"<<endl;
	  
	  k++;
	  
	  fStationList.push_back(*sta);
	  delete sta;
	  //cout<< "station = "<< k <<endl;
	}
      }
    
    fNorMax = NMIN + n * 1300;
    fNorMin = NMIN;  
    fEasMax = EMIN + n * 750 ;
    fEasMin = EMIN;

  }
  else{
    
    
    n=(Int_t)sqrt((float)size);
    
    // cout<<"Write a default array of "<<size<<" stations"<<endl;
    for(Int_t i=0;i<n;i++)
      {
	for(Int_t j=0;j<n;j++) 
	  {
	    Station * sta= new Station (k,NMIN+i*1609,EMIN+j*1609,1200);
	    stafile<<k<<"\t"<<NMIN+i*1609<<"\t"<<EMIN+j*1609<<"\t 1200"<<endl;
	    
	    k++;
	    
	    fStationList.push_back(*sta);
	    delete sta;
	    //cout<< "station = "<< k <<endl;
	  }
      }
    
    fNorMax = NMIN + n * 1609;
    fNorMin = NMIN;  
    fEasMax = EMIN + n * 1609 ;
    fEasMin = EMIN;
  }
  fNorMean = (fNorMax+fNorMin)/2.;
  fEasMean = (fEasMax+fEasMin)/2.;
  fNorWidth = fNorMax-fNorMin ;
  fEasWidth = fEasMax-fEasMin;

  stafile.close();

}

Array::Array(string arrayfilename)
{
  Int_t id;
  Double_t north,east;
  string name,endline;
  ifstream arrayfile;
  fNorMin=99999999;
  fNorMax=0;
  fEasMin=99999999;
  fEasMax=0;

  cout<<endl<<"Reading array in file = "<< arrayfilename<<endl<<endl;
  
  int nsta=0;
  
  arrayfile.open(arrayfilename.c_str(),ios::in);
  
  if(arrayfile==0)	
    {	
      cout << "++ error in opening file " << arrayfilename.c_str()  << endl;
      exit(1);
    }		
  
  cout<<"stations coordinates in Northing Easting mode"<<endl;
 
  while(!arrayfile.eof())
    {      
      arrayfile>>id>>east>>north;
      nsta++;
      cout<<id<<" id "<<north<<" entre "<<east<<" ini"<<endl;

      Station *sta =new Station(id,north,east);
      fStationList.push_back(*sta);
      delete sta;
      //cout<< "station = "<<id <<endl;
      
      if(north>fNorMax)fNorMax = north;
      if(north<fNorMin)fNorMin = north;  
      if(east>fEasMax)fEasMax = east;
      if(east<fEasMin)fEasMin = east; 
    }
    
  fNorMean = (fNorMax+fNorMin)/2. ;
  fEasMean = (fEasMax+fEasMin)/2. ;
  fNorWidth = fNorMax-fNorMin ;
  fEasWidth = fEasMax-fEasMin ;
  cout<<"array build with "<<nsta<<" stations"<<endl;
  cout<<fNorMin<<" "<<fNorMax<<" "<<fEasMin<<" "<<fEasMax<<endl;
}

Array::Array(Double_t theta,Double_t phi,Double_t E0,Double_t N0)
{
  const Int_t NBCIRCLES=99;
  const Double_t MINRADIUS=50.;
  const Double_t NBPOINTSINCIRCLE=12;
  const Double_t STEPINRADIUS=50;
  Double_t r;
  Double_t east;
  Double_t north; 
  Double_t altitude;
  fNorMin=99999999;
  fNorMax=0;
  fEasMin=99999999;
  fEasMax=0;
  
  Int_t nbstations;
 
  ofstream stafile;
  stafile.open("stararray.txt",ios::out);

  for(Int_t i=0; i< NBCIRCLES; i++)
    {      
      r= (Double_t)(MINRADIUS + i * STEPINRADIUS);
      for(Int_t j=0; j<  NBPOINTSINCIRCLE;j++)
	{
	  nbstations++;
	
	  altitude=0.;
	  //	Double_t angle=Double_t(j+.5)*2*M_PI/Double_t(NBPOINTSINCIRCLE);
	  Double_t angle=Double_t(j)*2.*PI/Double_t(NBPOINTSINCIRCLE);
	  Double_t xsh=r*cos(angle)/cos(theta);
	  Double_t ysh=r*sin(angle);
	  Int_t id=100*i+j+1;
	  //  east=E0+xsh*cos(theta)*cos(phi)-ysh*sin(phi);
	  // north=N0+ xsh*sin(phi)*cos(theta)+ysh*cos(phi); 
	
	  east=E0+xsh*cos(phi)-ysh*sin(phi)+altitude*tan(theta)*cos(phi);
	  north=N0+ xsh*sin(phi)+ysh*cos(phi)+altitude*tan(theta)*sin(phi); 
	             
	  Station * sta= new Station (id,north,east,ALTITUDE_SITE);
	  stafile<< id << " " << xsh << " "<< ysh <<" "<< east << " " << north << endl;
	       
	  fStationList.push_back(*sta);
	  delete sta;	  
	}
    }
  
  
  if(north>fNorMax)fNorMax = north;
  if(north<fNorMin)fNorMin = north;  
  if(east>fEasMax)fEasMax = east;
  if(east<fEasMin)fEasMin = east; 
  fNorMean = (fNorMax+fNorMin)/2. ;
  fEasMean = (fEasMax+fEasMin)/2. ;
  fNorWidth = fNorMax-fNorMin ;
  fEasWidth = fEasMax-fEasMin ;
	
  stafile.close();
    
}

void Array::Fill(Station* sta)  
{
  fStationList.push_back(* sta); 
}

Int_t Array::Size()
{
  Int_t size;
  size = fStationList.size();
  cout<<endl;
  cout<<"\tArray of " <<size << " stations"<<endl;
  return size;
}
