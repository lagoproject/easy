//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File Station.h 

 */
//----------------------------------------------------------------------

#ifndef STATION_H
#define STATION_H

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include "TObject.h" 
#include "Constants.h" 
#include "Particle.h" 

//----------------------------------------------------------------
/*
  class Station (parent class)
  A Station is defined by its coordinates.
*/
//----------------------------------------------------------------

  
class Station
{
 private:
  
 public:
  
  Int_t fId;
  string fName;
  Float_t fNorthing;
  Float_t fEasting;
  Float_t fAltitude; 
  Float_t fX_sf;
  Float_t fY_sf;
  Float_t fZ_sf;
  Float_t fR_sf;
  Float_t fAzim_sf;

  
  Station();
  Station(Int_t id, Float_t north,Float_t east,Float_t alt);   
  Station (Int_t id, Float_t north,Float_t east); //used to have an array of event cores
  virtual ~Station();

/*   void WriteName(string name); */
/*   void WriteShowerFrame(Double_t xstat,Double_t ystat,Double_t zstat); */
  // cartesian coordinates
  double X();
  double Y();
  double Z();
  // shower frame coordinates
  double X_sf();
  double Y_sf();
  double Z_sf();
  double R_sf();
  double Azim_sf();
  
  ClassDef(Station,1)
};


//----------------------------------------------------------------
/*
  class SampStation
  A SampStation contains all particles that have to be simulated
*/
//----------------------------------------------------------------


class SampStation : public Station
{
 private:

 public:
  
  Int_t fNmu;
  Int_t fNph;
  Int_t fNel;
  Int_t fNPartTot; 
  vector<Particle> fPartList;
   
  SampStation();
  SampStation(Int_t id, Float_t north,Float_t east,Float_t alt);   
  SampStation (Int_t id, Float_t north,Float_t east); //used to have an array of event cores
  virtual ~SampStation();
 
  void AddEntry(Particle *p);
 
  ClassDef(SampStation,1)
};


//----------------------------------------------------------------
/*
  class PMT
  A PMT.
*/
//----------------------------------------------------------------

class PMT
{
 private:

 public:

  Float_t fNpe;
  Float_t fNpe_direct;
  Float_t fNpe_mu;
  Float_t fNpe_ph;
  Float_t fNpe_el;
  Int_t fFirstPelTime;
  Int_t fFirstPelBin;
  Int_t fFirstADCBin;
  Double_t fSampFact;
  map<Int_t,Double_t> fTimeProfile;
  map<Int_t,Double_t> fTimeProfile_mu;
  map<Int_t,Double_t> fTimeProfile_em;
  map<Int_t,Double_t> fPMTSignal_hi;
  map<Int_t,Double_t> fPMTSignal_lo;
  Double_t fFEFilter_hi[MAXNUMBEROFTIMEBINS];
  Double_t fFEFilter_lo[MAXNUMBEROFTIMEBINS];
  short fADC[2][MAXNUMBEROFADCBINS];
  int fADC_ns[MAXNUMBEROFADCBINS];
  short fADC_mu[2][MAXNUMBEROFADCBINS];
  short fADC_em[2][MAXNUMBEROFADCBINS];

  PMT();
  virtual ~PMT();
  //  PMT(const PMT& p);

  void Init();
  void WritePE(Int_t idpart,Int_t time);
  void DirectLight(int npe_direct);
  void DoElecSim();
  void PMTAmplification(map<Int_t,Double_t>& spec ,map<Int_t,Double_t>& specAn,map<Int_t,Double_t>& specDy);
  void FEFilter(map<Int_t,Double_t>& spec1,Double_t* spec2);
  void FADCSampling(Double_t* spec1 ,short* spec2);
  void FADCSamplingSat(Double_t* spec1 ,short* spec2,int* spec3);
  

  ClassDef(PMT,1)

};


//----------------------------------------------------------------
/*
  class HitStation
  A HitStation contains all informations after simulation
*/
//----------------------------------------------------------------

class HitStation : public TObject,public SampStation
{
 private:

 public:

  Int_t fT2ToT;
  Int_t fT1Threshold;
  Int_t fT2Threshold;
  Int_t fTriggerBin;
  Int_t fFirstPelTime;
  Int_t fFirstPelBin;
  Double_t fTime0;
  Int_t fNtop;
  Int_t fNside;
  Float_t fNpe;
  Float_t fNpe_direct;
  Float_t fNpe_mu;
  Float_t fNpe_ph;
  Float_t fNpe_el;
  Double_t fSignalInVem;
  Double_t fSampFact;
  Int_t fNbOfBinsOverThres;
 vector<Int_t> fMuTimes;
  map<Int_t,Double_t> fTimeProfile;
  map<Int_t,Double_t> fTimeProfile_mu;
  map<Int_t,Double_t> fTimeProfile_em;
  short fADC[2][MAXNUMBEROFADCBINS];
  short fADC_ns[MAXNUMBEROFADCBINS];
  short fADC_mu[2][MAXNUMBEROFADCBINS];
  short fADC_em[2][MAXNUMBEROFADCBINS];
  Double_t fDoubleADC[2][MAXNUMBEROFADCBINS];
  Double_t fDoubleADC_mu[2][MAXNUMBEROFADCBINS];
  Double_t fDoubleADC_em[2][MAXNUMBEROFADCBINS];
  PMT fPMT[NPM];

  HitStation();
  HitStation(Int_t id, Float_t north,Float_t east,Float_t alt);   
  HitStation (Int_t id, Float_t north,Float_t east); //used to have an array of event cores
  ~HitStation();
  // HitStation(const HitStation& s);
  // HitStation& operator=(const HitStation& s);

  void WriteSamplingFactor(Double_t sampfact);
  void WritePE(Int_t idpart,Int_t time);
  void WriteADCBin(Int_t idpart,Int_t bin,Double_t bincontent);
  void WriteDoubleADCBin(Int_t idpart,Int_t bin,Double_t bincontent);
  void CountPart(Int_t idpart,Double_t wpart_det,Double_t time);
  void SumADCTraces();
  void DoDetSim(string simmode,double costheta_sh,double sintheta_sh,
		double cosphi_sh, double sinphi_sh,const SampStation &);
  void WriteShortADCTraces();
  void DoElecSim();
  void SetTime0();
  void SetTime0(double time0);
  double SignalInVem();
  void ReadParticleInTank(ParticleInTank* p);
  int DoTrigSim();
  int TestToT();
  int TestThreshold(int trigtype);
  ClassDef(HitStation,1)

};


#endif
