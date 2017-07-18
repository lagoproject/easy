using namespace std;

#include "Constants.h"
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

int gNbins;
double gNpesamples;
int gNshapezones;
int gNshape;
int gShapezone[MAXNUMBEROFADCBINS];
double gIntshape[MAXNUMBEROFADCBINS];

 
double MuonPulseShape(int ti)
{
  double value;
  if(ti<=0)value =0;
  else
    value = - MUONPULSE_RISINGNORM * exp(-(ti+8)/ MUONPULSE_RISINGTIME) +
      MUONPULSE_DECAYNORM * exp( -(ti+8)/ MUONPULSE_DECAYTIME );
  return(value);
}

void BuildIntegratedMuonPulse(string mode)
{
  double totaltime;
  
  
//  Sampling time related quantities.
  if(mode =="CALIB")
    totaltime= MUONPULSE_TOTALTIME_CALIB;
  else
    totaltime= MUONPULSE_TOTALTIME_SHOWER;
  
  gNbins = (int) (totaltime / SAMPLINGTIME);
  
  if (gNbins < 10) gNbins = 10;
  if (gNbins > MAXNUMBEROFADCBINS) {
    cerr<<"Too many bins implied by totaltime/samplingtime specification."<<endl;
    exit(0);
  }

  //     Calculating the number of bins needed to sample adequately
  //     a muon pulse.

  
  int ti;
  double shape[MAXNUMBEROFADCBINS];
 

  for(int i=0;i<MAXNUMBEROFADCBINS;i++){
    shape[i]=0;
    gIntshape[i]=0;
  }
 
  gNshape=0;
  shape[0] = 0;
  shape[1] = 0;
  double testvalue = -1;


  for(int i=0;i<MAXNUMBEROFADCBINS;i++)
    {
      gNshape++;
      shape[1] = MuonPulseShape(gNshape);
  
      if (shape[1] >= shape[0] )  
	{
	  shape[0] = shape[1];
	  testvalue = shape[0]/1001;
	 
	}
      else if (testvalue > shape[1]  )
	break;
     
    }


  
  gNshape = (int) (0.5 + (double) gNshape / SAMPLINGTIME);

  if (gNshape < 1) gNshape = 1;
  if (gNshape >MAXNUMBEROFADCBINS ) 
    {
      cerr<<   "Too many shape samples."<<endl;
      exit(0);
    }
  
  
  //     Evaluating normalized cumulative function, and number of samples
  // for each photoelectron.
  
  
  
  for(int i=0;i<gNshape;i++)
    {
      ti =(int)(((double)i-0.5)* SAMPLINGTIME);
      shape[i] = MuonPulseShape(ti);
      if(i >0)gIntshape[i] = shape[i]+gIntshape[i-1] ;
     
    }
  
  gNpesamples = 1;
  
  for(int i=0;i<gNshape;i++)
    {
      gIntshape[i] /= gIntshape[gNshape-1];
      if (gIntshape[i] < MUONPULSE_PROBLIMIT) gNpesamples = i;
    }
  
  cout<<"Muon pulse: Shaping with "<< gNshape<<" bins ("<< gNshape * SAMPLINGTIME<< " ns)"<<endl;
  cout<< "Muon pulse: Integrated 75% at bin "<< gNpesamples <<
    " ("<< gNpesamples * SAMPLINGTIME<< " ns)"<<endl;
  
  
  
  //     Evaluating zone index for (fast) indexed search of shape array.
  double sht,shi;
  
  gNshapezones = gNshape - 1;
  sht         = gNshapezones;

 
  gShapezone[0] = 1;
  int k;
  for(int i=0;i<gNshapezones;i++)
    {
      shi = i / sht;
      for(int j=0;j<gNshape;j++)
	{
	  k = j;
          if (gIntshape[j] >= shi) break;
	}
      gShapezone[i] = k;
      
    }
  
}
