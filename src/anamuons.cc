#include "TROOT.h"
#include "TRint.h"
#include"TFile.h"
#include "TRandom.h"
#include "Analyze.h"
using namespace std;

TFile* inrootfile,*outrootfile;
//string inrootfilename; 
//string outrootfilename; 
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------






//int Error; // needed by motif
extern void InitGui();//initializer for GUI needed for interactive interface
VoidFuncPtr_t Initfuncs[]={InitGui,0};

// initialize the Root system
//TROOT root("Rint","The SDS-ROOT interactive Interface",Initfuncs);
TROOT root("Rint","Analyze of muons");

int main(int argc,char** argv)
{
  
  //  create the interactive interface
  TRint *theApp=new TRint("ROOT-SDSIM",&argc,argv);
  theApp->Run();
  return 0;
}


