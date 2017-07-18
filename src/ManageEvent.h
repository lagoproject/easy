//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File ManageEvent.h 

 */
//----------------------------------------------------------------------

#ifndef MANAGEEVENT_H
#define MANAGEEVENT_H

#include <TROOT.h> 
#include "Event.h" 
#include "Station.h" 
#include "Constants.h" 

//----------------------------------------------------------------
/*
  Functions to manage a simulated event in mode SHOWER
*/
//----------------------------------------------------------------

Event *GenerateEvent(UInt_t);
void DoSampling(Double_t, Double_t, vector<SampStation>*);
void AddRandomParticles(vector<SampStation>*);
HitStation *GenerateHitStation(const SampStation &);

//----------------------------------------------------------------
/*
  Functions to manage a simulated event in mode CALIB
*/
//----------------------------------------------------------------


#endif
