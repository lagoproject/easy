#ifndef _ANALYZESTAR_H
#define _ANALYZESTAR_H


#include "TObject.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include "TBrowser.h"
#include "Event.h"
#include "Array.h"
#include "TBranch.h"
#include <vector>
#include <map>
#include "TChain.h"





class AnalyzeStar : public TObject{
 private:
  TList *fCanvasList;
  TList *fHisto1;
  TList *fHistTitle;
  TList *fHistAxis;
  TCanvas* fCanvas;
  TFile* fFile;
  TTree * fTree;
  TChain* fChain;
  Event* fEvent;
  TBranch* fBranch;
  Array* fArray;
  Array* fCoreArray;
  HitStation* fHitStation;
  vector<HitStation> fStList;
 
  
 public:
  AnalyzeStar();
  ~AnalyzeStar();
  void CompareProfiles();
  void ComputeNbotDistr(int energy, int tta);
  void DoLadybot(); 
  void DoLadybotToT();
  void DoLTP(char *energy);
  void DrawSaturatedTrace(int ipm );
  void DrawSaturatedTrace();
  void GetStation(Int_t sta);
  void GetStation(Int_t sta,Int_t gain);
  void LateralProfile();
  void PlotNbotDistr();
  void PrintS1000(char *energy);
  void PrintTree();
  void ReadFilesCDF(char * Energy,char * tta);
  void ReadFilesLTP();
  void ReadFile(char * file ,int opt);
   
   /*
 
  void DrawEvent();
  void Help();
  */
 private:
 

  ClassDef(AnalyzeStar,1)
    };
#endif
    







