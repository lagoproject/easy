
/*
  ShowerParam.cc

  implementation file for class ShowerParam

*/


#include "EasySim.h"
#include "ShowerParam.h"
#include <TROOT.h>
#include <TFile.h>
#ifndef CALIBONLY
#include "Shower_ROOT.h"
#include "ShowerGrnd_ROOT.h"
#include "PartGrnd_ROOT.h"
#include "Header_ROOT.h"
#include "HeaderAires_ROOT.h"
#endif
#include "Constants.h"


//----------------------------------------------------------------
/*
  class ShowerParam
  The input shower.
*/
//----------------------------------------------------------------

ShowerParam *gShParamP;

ShowerParam::ShowerParam()
{

}

/*
  Read the grdpcles files converted in the root format from Aires or corsika
  Rewrites the particles to a new format that includes a rotation angle that is
  needed when rotating a shower.
  Note : theConfig()->InShower is allocated here
*/
ShowerParam::ShowerParam(string file)
{
#ifndef CALIBONLY
  TFile *rootfile;
  TObject *obj;

  char * filename;
  Int_t l = file.length()+1;
  filename = new char[l];

  rootfile=new TFile(file.c_str(),"R");

  //  rootfile->ls();

  obj=(Shower_ROOT*)rootfile->Get("Header_ROOT");

//   if( obj!=0)
//     {
//       cerr << "\t Header_ROOT object found " << endl;
//     }

  //Reads shower characteristics
  Primary=(Int_t)((Header_ROOT*)obj)->GetPrimaryType();
  Energy=((Header_ROOT*)obj)->GetPrimaryEnGev();
  Theta=((Header_ROOT*)obj)->GetPrimaryZenAng();
  Phi=((Header_ROOT*)obj)->GetPrimaryAzAng();

  cout<<"\tSHOWER CHARACTERISTICS ==>"<<endl<<endl;
  cout<< "\t\t\tPrimary = "<< Primary<<endl<<"\t\t\tEnergy = " 
      << floor(Energy*1.E-9*100)/100. <<" EeV"<<endl
      << "\t\t\tTheta = "<<floor(Theta*RAD2DEG*100)/100.<<" °"<<endl<<"\t\t\tPhi= " 
      << floor(Phi*RAD2DEG*100)/100.<<" °"<<endl;

  ShowerProf_ROOT* profile=0;

  if( (profile=(ShowerProf_ROOT*)rootfile->Get("ShowerProf_ROOT"))!=NULL)
    {
      XMax=profile->GetChimax();
      X0=profile->GetChi0();
    }

  strcpy(filename,file.c_str());

  theConfig()->InShower=new Shower_ROOT(filename);
  NParticles=(Int_t)theConfig()->InShower->showerGrnd()->GetEntries();
  cerr<<endl;
  cout<<"\t==> "<<NParticles<<" particles to sample" <<endl<<endl;

  cout<<endl;for (int ii=0; ii<80; ii++) cout<<"-";cout<<endl;

  if(theConfig()->PhiRotation!=0.) Phi +=  theConfig()->PhiRotation * DEG2RAD;
  CosPhi = cos(Phi);
  SinPhi = sin(Phi);
  CosTheta = cos(Theta);
  SinTheta = sin(Theta);

  delete[] filename;
  delete rootfile;
#endif
}

ShowerParam::~ShowerParam()
{

}
