//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File ShowerParam.h 

 */
//----------------------------------------------------------------------

#ifndef SHOWERPARAM_H
#define SHOWERPARAM_H

#include <TROOT.h>
#include <string>

//----------------------------------------------------------------
/*
  class ShowerParam
  The input shower.
*/
//----------------------------------------------------------------

class ShowerParam
{
 private:

 public:

  Int_t Primary;
  Int_t NParticles;
  Double_t Theta;
  Double_t Phi;
  Double_t CosTheta;
  Double_t CosPhi;
  Double_t SinTheta;
  Double_t SinPhi;
  Double_t Energy;
  Double_t XMax;
  Double_t X0;
  Double_t XCore;
  Double_t YCore;
  Double_t EasCore;
  Double_t NorCore;
  
  ShowerParam();
  ShowerParam(string file);
  ~ShowerParam();
  

};
extern ShowerParam *gShParamP;


#endif
