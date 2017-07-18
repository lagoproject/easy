#include "Analyze.h"
#include <iostream>
using namespace std;

int main(int argc,char** argv) {
  if (argc!=2) {
    cerr << "Syntax: " << argv[0] << " file.root" << endl;
    cerr << "  Produces an ascii file with the following format:" << endl;
    cerr << "  ParticuleNumber FADCBin(0-100) PMId(0-n) ADCValue" << endl;
    return 1;
  }
  TFile * File = new TFile(argv[1],"READ"); 
  TTree* Tree = (TTree *) File->Get("TreeEvent");
  TBranch* Branch= Tree->GetBranch("event");
  Event* Event = 0;
  Tree->SetBranchAddress("event",&Event);

  HitStation* sta;
  vector<HitStation> StList;
  int NEntries=(int)Tree->GetEntries();
  cerr << NEntries << " particles to analyze" << endl;
  cout << "# Ascii dump version 0.2" << endl;
  cout << "# 2" << endl;
  cout << "# ParticuleNumber FADCBin(0-100) PMId(0-n) ADCValue" << endl;
  for(int i=0; i<NEntries; i++) {
    Tree->GetEvent(i);
    if(i%100==0) cerr << i << " particles done" << endl;
    StList = Event->fHitStationList;
    sta=&(StList[0]);
    for (int ipm=0;ipm<NPM;ipm++)
      for(int k=0; k<100; k++)
        cout << i << " " << k << " " << ipm << " " << sta->fPMT[ipm].fADC[0][k] << endl;
  }
  cerr << NEntries << " particles done" << endl;
  return 0;
}

