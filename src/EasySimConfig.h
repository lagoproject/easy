//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File EasySimConfig.h 

 */
//----------------------------------------------------------------------

#ifndef EASYSIMCONFIG_H
#define EASYSIMCONFIG_H

#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "Utils.h"
#ifndef CALIBONLY
#include "Shower_ROOT.h"
#include "ShowerGrnd_ROOT.h"
#include "PartGrnd_ROOT.h"
#include "Header_ROOT.h"
#include "HeaderAires_ROOT.h"
#endif

//----------------------------------------------------------------
/*
  class EasySimConfig
  The configuration parameters of a run.
*/
//----------------------------------------------------------------

class EasySimConfig
{
 private:
  static EasySimConfig *_instance;

 protected:
  EasySimConfig();

 public:

  string RunNumber;
  string Mode;
  string SimMode;
  string MuEmMode;
  string ElecMode;
  string PartMode;
  string PartAddMode;
  Int_t NEvents;
  Int_t PartCode;
  Double_t PartEnergy;
  Double_t PartTheta;
  Double_t EasAvrg;
  Double_t EasSpread;
  Double_t NorAvrg;
  Double_t NorSpread;
  Int_t NPartMultiple;
  string ArrayMode;
  string ArrayFileName;
  string EventMode;
  Double_t PhiRotation;
  Double_t RMin;
  Double_t RMax;
  Double_t DeltaR;
  Double_t DeltaAzim;
#ifndef CALIBONLY
  Shower_ROOT *InShower;  
#endif
  TFile *OutRootFile;
  TTree *OutTree;
  TTree *ArrayTree;
  TTree *CoreTree;
  TBranch *OutBranch;
  TBranch *ArrayBranch;
  TBranch *CoreBranch;
  string InRootFileName;
  string OutRootFileName;
  string AsciiFileName;
  Int_t NoFileName;
  Int_t TwoFileName;
  Int_t AsciiFile;
  Int_t MuEmFlag;
  Int_t Acceptance;

  static EasySimConfig* Instance();
  
};
EasySimConfig* theConfig();


#endif
