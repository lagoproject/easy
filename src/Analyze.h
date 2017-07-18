#ifndef _ANALYZE_H
#define _ANALYZE_H


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





class Analyze : public TObject{
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
  Analyze();
  ~Analyze();

  void DirectLight();
  void DoLadybot3D();       
  void DrawArray();
  void DrawAnodeDynode();
  void DrawBalance(int ipm);
  void DrawEvent(int nev);
  void DrawMuonDecay();
  void DrawMuons();
  void DrawMuons(int opt);
  void DrawMuons(int n1,int n2,int opt);
  void DrawMuonSpectrum();
  void DrawMuonsCharge(int ipm);
  void DrawMuonsFE();
  void DrawMuonTimes();
  void DrawOneFE(int num,int ipm);          //Draw sum of PMTS
  void DrawOneFE(int num,int ipm,int ichan);          //Draw sum of PMTS
  void DrawOnePM(int num,int ipm);          //Draw sum of PMTS
  void DrawOnePM(int num,int ipm,int ichan);          //Draw sum of PMTS
  void DrawOneProfile(int num);          //Draw sum of PMTS
  void DrawOneProfile(int num,int ipm);          //Draw sum of PMTS
  void DrawOneTrace(int num);          //Draw sum of PMTS
  void DrawOneTrace(int num,int ipm);
  void DrawOneTrace(int num,int ipm,int ichan);
  void DrawOneTraceInVEM(int num,int ipm,int ichan);
  void DrawPulseShape();
  void DrawTraces();                   //Draw sum of PMTS
  void DrawTraces(int ipm);
  void DrawTracesInVem();                   //Draw sum of PMTS
  void DrawTracesInVem(int ipm);
  void GetEvent(int event);
  void GetEvents();
  void GetEvents(int ista);
  void PlotDistances();
  void PlotPMTBalance();
  void PrintTree();
  void ReadFile(char * file );
  void TestAnodeDynode();
  void TestNewTOT();
 
  /*
 
  void DrawEvent();
  void Help();
  */
 private:
 

  ClassDef(Analyze,1)
    };
#endif
    







