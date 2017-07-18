//============================================================================ 
//
// Header-File : Constants.h 
// Description : This header file declares some constants used in the package
// 
//============================================================================

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <string>
//----------------------------------------------------------
// Global Variables:
//----------------------------------------------------------

//GENERAL CONSTANTS

const double EPSILON            = 1.e-6;
const double PI                 = 3.14159265358979323846;

const double DEG2RAD            = PI/180;
const double RAD2DEG            = 180/PI;
const double SQRT3              = 1.732050808;
const double CLIGHT             = 0.299792458;// Speed of light given in m / ns
const double GFERMI             = 1.16639E-5; // GeV^-2

//LAGO SITE CONSTANTS
const double ALTITUDE_SITE		= 865.; //(m) reference altitude

//TANK CONSTANTS
const double STATION_RADIUS     = 0.78;  //(m)
const double STATION_HEIGHT     = 1.54;  //(m)
const int  NPM                  = 1;
const double X_PM[NPM]          = {0.};
const double Y_PM[NPM]          = {0.};
const double RAD_PM             = .1477;//(m)
const double HC_PM              = .0776;//(m)
const double TOP_FACT           = 1.;  //white top

//SIMULATION CONSTANTS
const int NCLONES =  10;
const int NMAXPARTINTANK  =5000;

//BINNING OF PE PROFILES AND ADC TRACES

const int  MAXNUMBEROFTIMEBINS  = 19200/*10000*/; //bins of 1ns
const int  SAMPLINGTIME         = 25; // (ns)
const int  MAXNUMBEROFADCBINS   = MAXNUMBEROFTIMEBINS/SAMPLINGTIME;
const int  TIMEOFFSET           = 100 ; //time offset of 100 ns ( to save very early particles) 
 


const double MASS_E              = .000511;  //( in GeV)
const double MASS_MU             = .1057;   //( in GeV)
const double LIFE_MUON           =  2197.; //  ( in ns)
const double EMIN_E              = .00025; //
const double EMIN_MU             = .05 ;   //
const double DEDX_E              = .2 ;    //   in GeV/m
const double DEDX_MU             = .2;     //in GeV/m
const double CHERCONST           = .0458506;
const double DELTARAYCONST       = 0.008445; // in Gev/m







//CONSTANTS FOR FAST SIMULATION BASED ON MUONS SHAPE

const double MUONPULSE_RISINGTIME       = 13.31;
const double MUONPULSE_DECAYTIME        = 54.72;
const double MUONPULSE_RISINGNORM       = 3.77;
const double MUONPULSE_DECAYNORM        = 2.45;
const double MUONPULSE_TOTALTIME_CALIB  = 400.; 
const double MUONPULSE_TOTALTIME_SHOWER = 10000.;
const double MUONPULSE_PROBLIMIT        = 0.75;
const double VEMNORMALIZATION           = 13.6; //to normalize the fast sim to the vem


#endif // __CONSTANTS_H

//==End of File=================================================================
