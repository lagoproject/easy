//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File Event.h 

 */
//----------------------------------------------------------------------

#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include <string>
#include "TObject.h"
#include "Station.h"

//----------------------------------------------------------------
/*
  class Event (parent class)
  A simulated event.
*/
//----------------------------------------------------------------

class Event : public TObject
{
 private:

 public:

  Int_t fId;
  Int_t fPrimary;
  Double_t fEnergy;
  Double_t fTheta;
  Double_t fAzim;
  Double_t fNorCore;
  Double_t fEasCore;
  Double_t fPhiRot;
  Double_t fXCore;
  Double_t fYCore;
  Double_t fZCore;
 Double_t fXmax;
  Double_t fX0;
  Int_t fNombTank;
  Int_t fT3Algo;
  Int_t fT4Algo;
  vector<HitStation> fHitStationList;

  Event();
  Event(Int_t primary,Double_t energy,Double_t theta ,Double_t azim);
 
   Event(Int_t primary,Double_t energy,Double_t theta ,Double_t azim,Double_t xcore ,Double_t ycore,Double_t zcore);
 Event(Int_t primary,Double_t energy);
  Event(Int_t primary,Double_t energy,Double_t theta ,Double_t azim,
        Double_t eas_avrg,Double_t eas_spread, Double_t nor_avrg, Double_t nor_spread);
  ~Event();
  //  Event(const Event& e);

  void CountTanks();
  void AddTank(HitStation *hst);
  void SetDepth(Double_t xmax,Double_t x0);
  int DoT3Sim();
  int DoT4Sim();
  ClassDef(Event,1)
};


#endif
