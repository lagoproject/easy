//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File Particle.h 

 */
//----------------------------------------------------------------------

#ifndef PARTICLE_H
#define PARTICLE_H

using namespace std;

#include "TObject.h"
#include <vector>
#include <string>

//----------------------------------------------------------------
/*
  class Particle (parent class)
  A particle at the ground level.
  This particle has be produced by an Air Shower simulation code.
*/
//----------------------------------------------------------------

class Particle : public TObject 
{
 private:

 public:

  Int_t fId;            /* Aires code */
  Float_t fUX;
  Float_t fUY;
  Float_t fUZ;
  Float_t fT;           /* T: time of arrival */
  Float_t fT_plane;     /* T: time of  */
  Double_t fE;           /* En: energy in Gev */
  //  Float_t fE;           /* En: energy in Gev */
  Float_t fWtop;
  Float_t fWside;

  Particle();
  Particle(Int_t idpart,Int_t ntop,Int_t nside,Float_t time_plane,
	   Float_t timepart,Float_t energy,Float_t cx,Float_t cy, Float_t cz );
  ~Particle();

  ClassDef(Particle, 1)

};

//----------------------------------------------------------------
/*
  class ParticleInTank

*/
//----------------------------------------------------------------

class ParticleInTank : public Particle
{
 private:
 
  /* step parameters for charged particles */
  Double_t Beta;//!
  Double_t Range;//!
  Double_t DxStep;//!
  Int_t IStep; //!

  /* step parameters for gammas */
  Double_t GLength; //!
  Double_t MeanIntLength; //!
  Double_t IntLength; //!
  Double_t ProbPair; 

 public:
  Int_t fISta;
  Int_t fMotherId;
  Double_t fX;
  Double_t fY;		/*X,Y : position in m */
  Double_t fZ;  
  Double_t fW;		/* W : weight */

  Int_t fNCherPhot;     /* Number of Cherenkov photons induced */
  vector<Double_t> fTCherPhot;//!
  vector<Double_t> fNReflexions;//!
  vector<Int_t> fIpm;         //!
  vector<Double_t> fBinContent;//!
 
  ParticleInTank();
  ParticleInTank (Int_t ist,Int_t id,Double_t x,Double_t y,Double_t z, 
		  Double_t cx, Double_t cy, Double_t cz,Double_t tim,Double_t en,Double_t weight);
 ParticleInTank (Int_t ist,Int_t id,Int_t motherid,Double_t x,Double_t y,Double_t z, 
		  Double_t cx, Double_t cy, Double_t cz,Double_t tim,Double_t en,Double_t weight);  virtual ~ParticleInTank();
  
 
  Double_t GeomLength();
  void DoDetSim(string simmode);
  void GammaInt( string simmode);
  void CherPhot(Double_t beta,Double_t range);
  Int_t Follow(Double_t x,Double_t y,Double_t z,Double_t cosx,Double_t cosy,
	       Double_t cosz,Double_t prb_refdif,Double_t prb_spec,Double_t abslen);
  void FastSim(Int_t npe,Double_t wgt);
  Double_t  GenerateDeltaRays(Double_t beta, Double_t xstep);
  void MultipleScattering(Double_t beta,Double_t xstep);
  Double_t  Bremstrahlung(Int_t istep, Double_t dxstep,Double_t beta);
  int MuonDecay();
  void PairProduction();
  void ComptonScattering();

  void DoDetSim();
  void ComputeStepParameters();
  int DoProcesses();
  void Propagate();
  int DoNextStep();


  ClassDef(ParticleInTank, 1)

};
 

//----------------------------------------------------------------
/*
  class OutgoingParticle

*/
//----------------------------------------------------------------

/* class OutgoingParticle : public Particle */
/* { */
/*  private: */

/*  public: */

/*   Int_t fISta; */
/*   Double_t fX; */
/*   Double_t fY;           */
/*   Double_t fZ; */
/*   Double_t fW;           */

/*   OutgoingParticle(); */
/*   OutgoingParticle (Int_t ist,Int_t id,Double_t x,Double_t y,Double_t z,  */
/* 		    Double_t cx, Double_t cy, Double_t cz,Double_t tim,Double_t en, */
/* 		    Double_t weight); */
/*   ~OutgoingParticle(); */


/*   ClassDef(OutgoingParticle,1) */

/* }; */

#endif
