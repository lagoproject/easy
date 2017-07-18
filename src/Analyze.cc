#include "Analyze.h"
#include "Constants.h"
#include "Trigger.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <string.h>
#include <map>
#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TText.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TArrayD.h"
#include "TF1.h"
#include "TF2.h"
#include "TBRIK.h"
#include "TTUBE.h"
#include "TTUBS.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TNode.h"
#include "TRotMatrix.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraph2D.h"

using namespace std;

ClassImp(Analyze)


  Analyze::Analyze()
{
  fCanvasList=new TList();
  fHisto1= new TList();
  fHistTitle=new TList();
  fHistAxis=new TList(); 
 
  fTree = 0;
  fEvent =0;
  fBranch =0 ;
  fArray=0;
  fCoreArray=0;
  fFile=0;

 }

Analyze::~Analyze()
{
  fHistAxis->Delete();
  fHistTitle->Delete();
  fHisto1->Delete();
  fCanvasList->Delete();
  
  
  delete fHistAxis;
  delete fHistTitle;
  delete fHisto1;
  delete fCanvasList;
}

void Analyze::DirectLight()
{ 
  double npe,npedirect,npesum[NPM],npesumdirect[NPM];
  //  double nref;
  int nmu;
  for(int i=0;i<NPM;i++){
    npesum[i]=0;
    npesumdirect[i]=0;
    
  }
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);

  HitStation* sta;
  TH1F* hnpe[NPM];
  TH1F* hnpedir[NPM];
  for(int ipm=0;ipm<NPM;ipm++){
  hnpe[ipm] =  new TH1F("NPE","NPE",200,0.,200);
  hnpedir[ipm] =  new TH1F("NPEDIR","NEPDIR",200,0.,200);

  }
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;
  
  
  
  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      sta=&(fStList[0]);
      
      // work with adcsignal  
      for(int ipm=0;ipm<NPM;ipm++){
	npe =sta->fPMT[ipm].fNpe;
	npedirect =sta->fPMT[ipm].fNpe_direct;
	npesum[ipm]+=npe;
	npesumdirect[ipm]+=npedirect;
	if(npe>0)hnpe[ipm]->Fill(npe);
	if(npedirect>0)hnpedir[ipm]->Fill(npedirect);
	
      }
      
      
      
    }
  fCanvasList->Add(new TCanvas("npe"," npe",600,800));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,3);
  for(int ipm=0;ipm<NPM;ipm++){
    fCanvas->cd(2*ipm+1);
    hnpe[ipm]->Draw("hist");
    fCanvas->cd(2*ipm+2);
    hnpedir[ipm]->Draw("hist");
    
  }
  
  for(int ipm=0;ipm<NPM;ipm++)
    cout<<"pm"<<ipm+1<<" "<<100*(double)(npesumdirect[ipm])/(double)(npesum[ipm])<<" % direct light"<<endl;
}


  void Analyze::DoLadybot3D()
{	
  
  TH2F* Ladybot3DPM1 = new TH2F("Ladybot3DPM1","Ladybot3DPM1",80,0.,4000.,120,0.,119.);
  //TH2F* Ladybot3DPM2 = new TH2F("Ladybot3DPM2","Ladybot3DPM2",80,0.,4000.,120,0.,119.);
  //TH2F* Ladybot3DPM3 = new TH2F("Ladybot3DPM3","Ladybot3DPM3",80,0.,4000.,120,0.,119.);
 
  
  TCanvas *CanLadybot1= new TCanvas("Ladybot3D","Ladybot3D PM1",600,800);
  //TCanvas *CanLadybot2= new TCanvas("Ladybot3D","Ladybot3D PM2",600,800);
  //TCanvas *CanLadybot3= new TCanvas("Ladybot3D","Ladybot3D PM3",600,800);
  
  CanLadybot1->SetFillColor(0);
  CanLadybot1->SetBorderMode(0);

  vector<HitStation>::iterator sta;  
  
  for (int i=0;i< fChain->GetEntries() ; i++) 
    {
      fChain->GetEntry(i);
      fStList = fEvent->fHitStationList;
      
      for(sta=fStList.begin();  sta!= fStList.end(); sta++)
	{
	  
	  Ladybot3DPM1->Fill(sta->fR_sf,sta->fNbOfBinsOverThres);
	  //Ladybot3DPM2->Fill(sta->fR_sf,sta->fPMT[1].fNbOfBinsOverThres,1.);
	  //Ladybot3DPM3->Fill(sta->fR_sf,sta->fPMT[2].fNbOfBinsOverThres,1.);
	  
	  //Ladybot3DPM2->Fill(sta->fR_sf,0.,1.);
	  //Ladybot3DPM3->Fill(sta->fR_sf,0.,1.);
	  
	}
    }
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  CanLadybot1->cd();
  Ladybot3DPM1->GetXaxis()->SetTitle("Distance to axis in m");
  Ladybot3DPM1->GetYaxis()->SetTitle("Nb of bins over threshold");
  Ladybot3DPM1->GetZaxis()->SetTitle("Count");
  Ladybot3DPM1->Draw("surf");
  
  TCanvas *Canproject= new TCanvas("Projection","Projection",600,800);
  Canproject->cd();
  TH1D* projection=Ladybot3DPM1->ProjectionY("NbofBinsOverThres",20,40,""); 
  projection->Draw();
  //  CanLadybot2->cd();
  //Ladybot3DPM2->Draw("lego");

  
  //CanLadybot3->cd();
  //Ladybot3DPM3->Draw("lego");
   
  return;
}


void Analyze::DrawAnodeDynode()
{
  
  // short adcsum[2][MAXNUMBEROFADCBINS];
 
  
  HitStation* sta;
  int intdyn,an,dyn;
  int nbsta;
  TH2F* histo;
  TGraph* gr;
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  int nev=(int)fTree->GetEntries();
  
  cout<<nev<<" events to use"<<endl;
  
  histo = new TH2F("andyn","anode dynode",400,0,400,10000,0,10000);
  /*for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  int antab[10000];
  int dyntab[10000];
  
  for(int i=0;i<10000;i++){
    antab[i]=0;
    dyntab[i]=0;
  }
  
  for(int ne=0; ne<nev; ne++)
    {
      if(ne==0)cout<<"premiere iteration"<<endl;
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      nbsta= fStList.size();
      //  cout<<nbsta<<endl;
      for(int ista=0;ista<nbsta;ista++){
	sta=&(fStList[ista]);
	
	
	if(sta->fNpe>0) 
	  {
	    // adcsum =sta->fPMT[0].fADC;
	    intdyn=0;
	    int iflag=0;
	    int first=0;
	    int flagsat=0;
	    
	    for(int k=0; k<100; k++){
	      if(sta->fPMT[0].fADC[0][k]!=0){
		intdyn+=sta->fPMT[0].fADC[0][k] ;
		if(iflag==0 && sta->fPMT[0].fADC[1][k]>=2)
		  {
		    first=k;
		    iflag=1;}
		if (sta->fPMT[0].fADC[0][k]==1023)flagsat=1;
	      }
	    }
	    
	    if(flagsat==1)cout<<"saturation "<<ne<<" "<<sta->fId<<endl;
	    //if(sta->fNpe>2000) 
	    if(intdyn>0 && flagsat==0)
	      {
		an=0;
		dyn=0;
		//cout<<"stations used for anode/dynode"<< sta->fId<<endl;
		for(int k=first; k<20+first; k++){
		  dyn+=sta->fPMT[0].fADC[0][k] ;
		  an+=sta->fPMT[0].fADC[1][k] ;
		}
		//cout<<sta->fId<<" "<<an<<" "<< dyn<<endl;
		if(an==0)cout<<"pb in the concept "<< sta->fId<< " "<<first<<endl;
		else
		  {
		    histo->Fill(an,dyn);
		    antab[ne]=an;
		    dyntab[ne]=dyn;
		    // cout<<antab[ne]<<" "<<dyntab[ne]<<endl;
		  }
	      }
	    
	  }//end of conditions in npe
      }//end of loop on stations
    }//end of loop on event
  gr = new TGraph(10000,antab,dyntab);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(21);
  gr->SetTitle("anode dynode");
  gr->SetMarkerSize(1);
  gr->Draw("AP");   
  
}

void Analyze::DrawArray()
{
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);
  TH2F* harray;
  TH2F* hcores;
  Station* sta;


  fTree = (TTree *) fFile->Get("TreeArray");      // Get Tree
  fTree->Print(); 
  fArray=0;
  
  fTree->SetBranchAddress("array",&fArray);
  
  cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
  
  fTree->GetEntry(0);
  double xmin = (fArray->fEasMin)-1500;
  double xmax = (fArray->fEasMax);
  double ymin = (fArray->fNorMin)-1500;
  double ymax = (fArray->fNorMax);
  
  cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<endl;
  harray =  new TH2F("Array","Array",1000,xmin,xmax,1000,ymin,ymax);
  hcores =  new TH2F("Cores","Cores",1000,xmin,xmax,1000,ymin,ymax);
  int nomsta= fArray->fStationList.size();
  cout<<nomsta<<endl;
  for (int i= 0; i<nomsta;i++)
    {
      sta= &(fArray->fStationList[i]);
     
      harray->Fill(sta->fEasting,sta->fNorthing);
    }
  

 
  
  fTree = (TTree *) fFile->Get("TreeCore");      // Get Tree
  
  fTree->Print();
  fCoreArray=0; 
  fTree->SetBranchAddress("corearray",&fCoreArray);
  
  cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
  
  fTree->GetEvent(0);
  int nomcore= fCoreArray->fStationList.size();

  TText *stid[nomcore];
  for (int i= 0; i<nomcore;i++)
    {
      sta=&(fCoreArray->fStationList[i]);
      hcores->Fill(sta->fEasting,sta->fNorthing);
      char* name = new char[3];
      sprintf(name,"%3d",sta->fId);
      stid[i]=new TText(sta->fEasting,sta->fNorthing,name);
      delete name;
    } 
 
  fCanvasList->Add(new TCanvas("array"," array ",600,600));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  harray->SetMarkerStyle(20);
  harray->SetMarkerColor(kBlue);
  harray->Draw();
  hcores->SetMarkerStyle(21);
  hcores->SetMarkerColor(kRed);
  hcores->Draw("same");
  for(int i=0;i<nomcore;i++){
    stid[i]->Draw();
  }

}
void Analyze::DrawBalance(int opt){
 
  int nmu,npts=0,npesum;   
  int sigtot,sig[NPM];
  HitStation* sta;   	
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;
  float x[nmu],y[nmu],z[nmu],bal1[nmu],bal2[nmu],bal3[nmu],sigt[nmu];
  float balmoy1,balmoy2,balmoy3;
  TH1F * h1= new TH1F("bal1","bal1",100,0,1);
  TH1F * h2= new TH1F("bal2","bal2",100,0,1);
  TH1F * h3= new TH1F("bal3","bal3",100,0,1);
  TH1F * h= new TH1F("sigtot","sigtot",5000,0,25000);
  balmoy1=0;
  balmoy2=0;
  balmoy3=0;
  for(int ne=0; ne<1000; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;     
      sta=&(fStList[0]);
  
     // work with adcsignal        
      sigtot=0;      
      for(int k=0; k<400; k++){
	if(sta->fADC[0][k]!=0)
	  sigtot+=sta->fADC[0][k];	  
      }	 
      for (int i=0;i<NPM;i++)
	{
	  npesum=0;
	  for(int k=0; k<400; k++)
	    if(sta->fPMT[i].fADC[0][k]!=0){	     
	      npesum+=sta->fPMT[i].fADC[0][k] ;
	    } 
	sig[i]=npesum;  
	}   

      //if(sigtot>0 && fEvent->fZCore>=1.1987){
	  if(sigtot>0){
	x[npts]=fEvent->fXCore;
	y[npts]=fEvent->fYCore;
	//	z[npts]=fEvent->fZCore;
	sigt[npts]=(float)sigtot;
	bal1[npts]=(float)sig[0];///(float)sigtot;
	//bal2[npts]=(float)sig[1]/(float)sigtot;
	//bal3[npts]=(float)sig[2]/(float)sigtot;
	//balmoy1+=(float)sig[0]/(float)sigtot;
	//balmoy2+=(float)sig[1]/(float)sigtot;
	//balmoy3+=(float)sig[2]/(float)sigtot;
	h1->Fill(bal1[npts]);
	//h2->Fill(bal2[npts]);
	//h3->Fill(bal3[npts]);
	//h->Fill(sigtot);
	npts++;
      }
    }
  //balmoy1/=(float)npts;
  //balmoy2/=(float)npts;
  //balmoy3/=(float)npts;
  //cout<<balmoy1<<" "<<balmoy2<<" "<<balmoy3<<endl;
 //  TGraph2D *gr1= new TGraph2D(npts-1,x,y,bal1);
  //TGraph2D *gr2= new TGraph2D(npts-1,x,y,bal2);
  //TGraph2D *gr3= new TGraph2D(npts-1,x,y,bal3);
  
  if(opt==1){
  //   TGraph2D *gr1= new TGraph2D(npts-1,d,z,bal1);
//     TGraph2D *gr2= new TGraph2D(npts-1,d,z,bal2);
//     TGraph2D *gr3= new TGraph2D(npts-1,d,z,bal3);
    TGraph2D *gr= new TGraph2D(npts-1,x,y,sigt);
    
    fCanvasList->Add(new TCanvas("balance_PM1","balance_PM1  ",400,400));
    fCanvas=(TCanvas *) fCanvasList->Last();  
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(1);
    gr->GetXaxis()->SetTitle(" X core (m)");
    gr->GetYaxis()->SetTitle(" Y core (m)");
    gr->GetZaxis()->SetTitle(" PMT balance");
    gr->Draw("surf2z");
    /*
    fCanvasList->Add(new TCanvas("balance_PM2"," balance_PM2 ",400,400));
    fCanvas=(TCanvas *) fCanvasList->Last();  
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(1);
    gr2->GetXaxis()->SetTitle(" X core (m)");
    gr2->GetYaxis()->SetTitle(" Y core (m)");
    gr2->GetZaxis()->SetTitle(" PMT balance"); 
    gr2->Draw("surf2z");  
    
    fCanvasList->Add(new TCanvas("balance_PM3"," balance_PM3 ",400,400));
    fCanvas=(TCanvas *) fCanvasList->Last();  
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(1);
    gr3->GetXaxis()->SetTitle(" X core (m)");
    gr3->GetYaxis()->SetTitle(" Y core (m)");
    gr3->GetZaxis()->SetTitle(" PMT balance");
    gr3->Draw("surf2z");
   
    fCanvasList->Add(new TCanvas("totalsignal"," totalsignal ",400,400));
    fCanvas=(TCanvas *) fCanvasList->Last();  
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(1);
    gr->GetXaxis()->SetTitle(" X core (m)");
    gr->GetYaxis()->SetTitle(" Y core (m)");
    gr->GetZaxis()->SetTitle(" total signal");    
//gr->Projet("xy");
	gr->Draw("surf2z");*/
 }

  TCanvas * c1= new TCanvas("balances","balances  ",400,600);
  c1->Divide(1,3);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  h3->Draw();
  TCanvas * c2= new TCanvas("sigtot","sigtot",400,400);
  h->Draw();
 }

void Analyze::DrawEvent(int nev)
{
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);
  gStyle->SetPalette(1);
 TH2F* hcores;
  
  HitStation* sta;

  vector<HitStation> stlist;
  double xmin,ymin,xmax,ymax;

  xmin=99999999;
  ymin=99999999;
  xmax =0;
  ymax =0;

  fTree = (TTree *) fFile->Get("TreeEvent");      // Get Tree
  fTree->Print(); 
  fArray=0;
  
  fTree->SetBranchAddress("event",&fEvent);
  
  cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
  
  fTree->GetEntry(nev);
  stlist=fEvent->fHitStationList;
  int nombtank=stlist.size();
  TText *stid[nombtank];
  double latestime=9999999;
  cout<<nombtank<<endl;
  for(int i=0;i<nombtank;i++){
    sta= &(stlist[i]);
    if(sta->fNorthing > ymax) ymax = sta->fNorthing;
    if(sta->fNorthing < ymin) ymin = sta->fNorthing; 
    if(sta->fEasting > xmax) xmax = sta->fEasting; 
    if(sta->fEasting < xmin) xmin = sta->fEasting;
    char* name = new char[3];
    sprintf(name,"%3d",sta->fId);
    stid[i]=new TText(sta->fEasting+200,sta->fNorthing,name);
    stid[i]->SetTextSize(0.02);
    delete name;
    if(sta->fTime0<latestime)latestime=sta->fTime0;
   
  }
  xmin-= 1500;  
  ymin-= 1300;
  xmax+= 1500;
  ymax+= 1300;
  
  cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<endl;
  hcores =  new TH2F("Cores","Cores",1000,xmin,xmax,1000,ymin,ymax);
  
  
 
    hcores->Fill(fEvent->fEasCore,fEvent->fNorCore);
    
  
  fCanvasList->Add(new TCanvas("array"," array ",600,600));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
TStyle palette;
  palette.SetPalette(1);
  palette.Draw();
  hcores->SetMarkerStyle(21);
  hcores->SetMarkerColor(kRed);
  hcores->Draw();
  for(int i=0;i<nombtank;i++){
    stid[i]->Draw();
    sta= &(stlist[i]);
    double size= 0.5 +log(sta->fNpe)/5.;
    cout<<size<<endl;
    double time=sta->fTime0;
    TMarker *pt = new TMarker (sta->fEasting, sta->fNorthing, 20);
   //  pt->SetMarkerStyle(20);
//     pt->SetMarkerSize(size);
//     pt->SetMarkerColor(50 - (int)((time-latestime)/latestime*25));
//     // 	cout<<time-latestime<<endl;
//     pt->Draw("z");
       if(sta->fT2ToT){
      pt->SetMarkerStyle(20);
      pt->SetMarkerSize(size);
      pt->SetMarkerColor(50 - (int)((time-latestime)/latestime*25));
	cout<<time-latestime<<endl;
      pt->Draw("z");
    }
    else 
      if(sta->fT1Threshold){
      pt->SetMarkerStyle(20);
      pt->SetMarkerSize(1);
      pt->SetMarkerColor(5);
      pt->Draw();
      }
      else{
      pt->SetMarkerStyle(4);
      pt->SetMarkerSize(1);
      pt->SetMarkerColor(1);
      pt->Draw();
      }
  }
  
  
  


}

void Analyze::DrawMuonDecay()
{
  
  //short adcsum[2][MAXNUMBEROFADCBINS];
  map<int,double> trace;
  map<int,double> pmtrace;
  map<int,double> pmtracetot;

  int nmu;
  double npesum,npemax;
  HitStation* sta;
 
  TH1F*hadcch[NPM];
  TH1F*hadcsumch;
  TH1F*hpech[NPM];
  TH1F*hpesumch;
  TH1F*hpmch[NPM];
  TH1F*hpmsumch;

  TH1F*hadcpk[NPM];
  TH1F*hadcsumpk;
  TH1F*hpepk[NPM];
  TH1F*hpesumpk;
  TH1F*hpmpk[NPM];
  TH1F*hpmsumpk;

  TH1F*hadcatopk[NPM];
  TH1F*hadcsumatopk;
  TH1F*hpeatopk[NPM];
  TH1F*hpesumatopk;
  TH1F*hpmatopk[NPM];
  TH1F*hpmsumatopk;
  
  TProfile *profadcsum;
  TProfile *profpmsum;
  TProfile *profpesum;
  TProfile *profadc[NPM];
  TProfile *profpm[NPM];
  TProfile *profpe[NPM];
  
 
  gStyle->SetOptStat(1001111); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);
 
  hadcsumch =  new TH1F("ADCsumcharge","ADCsumcharge",2500,0.,5000);
  hadcsumpk =  new TH1F("adcsumpeak","adcsumpeak",200,0.,200);
  hadcsumatopk =  new TH1F("adcsumatopk","adcsumatopk",100,0.,10);
  profadcsum = new TProfile("profadcsum","profadcsum",400,0,400,0,450);
  hpesumch =  new TH1F("pesumcharge","pesumcharge",1500,0.,3000);
  hpesumpk =  new TH1F("pesumpeak","pesumpeak",100,0.,100);
  hpesumatopk =  new TH1F("pesumatopk","pesumatopk",100,0.,50);
  profpesum = new TProfile("profpesum","profpesum",8000,0,8000,0,50);
  hpmsumch =  new TH1F("pmsumcharge","pmsumcharge",2000,0.,500);
  hpmsumpk =  new TH1F("pmsumpeak","pmsumpeak",200,0.,1);
  hpmsumatopk =  new TH1F("pmsumatopk","pmsumatopk",200,0.,200);
  profpmsum = new TProfile("profpmsum","profpmsum",8000,0,8000,0,50);

  
   
  
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      char * hname1 = new char[15];
      sprintf(hname1,"ADCcharge_%1d",ipm+1);
      hadcch[ipm] =  new TH1F(hname1,hname1,100,0.,600);
      delete[] hname1;
      char* hname2 = new char[15];
      sprintf(hname2,"ADCpeak_%1d",ipm+1);
      hadcpk[ipm] =  new TH1F(hname2,hname2,100,0.,100);
      delete[] hname2;
     char* hname3 = new char[15];
      sprintf(hname3,"ADCatopeak_%1d",ipm+1);
      cout<<hname3<<endl;
      hadcatopk[ipm] =  new TH1F(hname3,hname3,100,0.,10);
       delete[] hname3;
      char* hname4 = new char[15];
      sprintf(hname4,"profADC_%1d",ipm+1);
      profadc[ipm] =  new TProfile(hname4,hname4,400,0.,400,0,50);
      delete[] hname4;
      char* hname5 = new char[15];
      sprintf(hname5,"PEcharge_%1d",ipm+1);
      hpech[ipm] =  new TH1F(hname5,hname5,200,0.,200);
      delete[] hname5;
      char*hname6 = new char[15];
      sprintf(hname6,"PEpeak_%1d",ipm+1);
      hpepk[ipm] =  new TH1F(hname6,hname6,50,0.,50);
      delete[] hname6;
      char*hname7 = new char[15];
      sprintf(hname7,"PEatopeak_%1d",ipm+1);
      hpeatopk[ipm] =  new TH1F(hname7,hname7,100,0.,50);
      delete[] hname7;
      char* hname8 = new char[15];
      sprintf(hname8,"profPE_%1d",ipm+1);
      profpe[ipm] =  new TProfile(hname8,hname8,8000,0.,8000,0,50);
      delete[] hname8;
      char*hname9 = new char[15];
      sprintf(hname9,"pmcharge_%1d",ipm+1);
      hpmch[ipm] =  new TH1F(hname9,hname9,100,0.,20);    
      delete[] hname9;
      char* hname10 = new char[15];
      sprintf(hname10,"pmpeak_%1d",ipm+1);
      hpmpk[ipm] =  new TH1F(hname10,hname10,200,0.,1);
      delete[] hname10;
      char* hname11 = new char[15];
      sprintf(hname11,"pmatopeak_%1d",ipm+1);
      hpmatopk[ipm] =  new TH1F(hname11,hname11,400,0.,200);
      delete[] hname11;
       char* hname12 = new char[15];
      sprintf(hname12,"profpm_%1d",ipm+1);
      profpm[ipm] =  new TProfile(hname12,hname12,8000,0.,8000,0,50);
      delete[] hname12;

    }

  

  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;
  
  
  
  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      sta=&(fStList[0]);
      
      // work with adcsignal  
      npesum=0;
      npemax=0;
      //   adcsum =sta->fADC;
       for(int k=0; k<400; k++){
	if(sta->fADC[0][k]!=0){
	  profadcsum->Fill(k,sta->fADC[0][k]);
	  npesum+=sta->fADC[0][k] ;
	  if(sta->fADC[0][k]>npemax) npemax=sta->fADC[0][k];
	  
	} 
	else
	  profadcsum->Fill(k,0);
      }// end of loop on k bins
      hadcsumch->Fill(npesum);
      hadcsumpk->Fill(npemax);
      if(npemax!=0)hadcsumatopk->Fill(npesum/npemax);
      
      
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  //  adcsum =sta->fPMT[ipm].fADC;
	  
	  // new loop on timebins
	  for(int k=0; k<400; k++){
	    if(sta->fPMT[ipm].fADC[0][k]!=0){
	      profadc[ipm]->Fill(k,sta->fPMT[ipm].fADC[0][k]);
	      npesum+=sta->fPMT[ipm].fADC[0][k] ;
	      if(sta->fPMT[ipm].fADC[0][k]>npemax) npemax=sta->fPMT[ipm].fADC[0][k];
	    } 
	    else
	      profadc[ipm]->Fill(k,0);
	  }// end of loop on k bins
	  hadcch[ipm]->Fill(npesum);
	  hadcpk[ipm]->Fill(npemax);
	  if(npemax!=0)hadcatopk[ipm]->Fill(npesum/npemax);
	  
	}
      
      npesum=0;
      npemax=0;
      
      trace=sta->fTimeProfile;
      
      for(int k=0; k<8000; k++){
	if(trace.count(k)){
	  npesum+=trace[k];
	  profpesum->Fill(k,trace[k]);
	  if(trace[k] > npemax) npemax=trace[k];
	}
	else
	  profpesum->Fill(k,0);
      }
      hpesumch->Fill(npesum);
      hpesumpk->Fill(npemax);
      if(npemax!=0)hpesumatopk->Fill(npesum/npemax);
      
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  
	  trace =sta->fPMT[ipm].fTimeProfile;
	  
	  for(int k=0; k<8000; k++){
	    if(trace.count(k))
	      {
	    
		npesum+=trace[k];
		profpe[ipm]->Fill(k,trace[k]);
		if(trace[k] > npemax) npemax=trace[k];
	      }
	    else
	      profpe[ipm]->Fill(k,0);
	      
	    }
	  hpech[ipm]->Fill(npesum);
	  hpepk[ipm]->Fill(npemax);
	  if(npemax!=0)hpeatopk[ipm]->Fill(npesum/npemax);
	  
	  
	}
      
      
      npesum=0;
      npemax=0;
      
      pmtracetot.clear();
      
      
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  
	  pmtrace =sta->fPMT[ipm].fPMTSignal_hi;
	  
	  for(int k=0; k<8000; k++){
	    if(pmtrace.count(k))
	      {
	      npesum+=pmtrace[k];
	      profpm[ipm]->Fill(k,pmtrace[k]);
	      if(pmtrace[k] > npemax) npemax=pmtrace[k];
	      if( pmtracetot.count(k))
		pmtracetot[k]+=pmtrace[k];
	      else
		pmtracetot[k]+=pmtrace[k];
	      if(pmtrace[k] > npemax) npemax=pmtrace[k];
	      
	    }
	    else
	      {
		profpm[ipm]->Fill(k,0);
	      }
	  }
	  hpmch[ipm]->Fill(npesum);
	  hpmpk[ipm]->Fill(npemax);
	  if(npemax!=0)hpmatopk[ipm]->Fill(npesum/npemax);
	    
	  
	}
      
      npesum=0;
      npemax=0;
      
      for(int k=0; k<8000; k++){
	if(pmtracetot.count(k))
	  {
	    npesum+=pmtracetot[k];
	    profpmsum->Fill(k,pmtracetot[k]);
	    if(pmtracetot[k] > npemax) npemax=pmtracetot[k];
	    
	    
	  }
	else
	  {
	    profpmsum->Fill(k,0);
	  }
      
      }
      for(map<int,double>::const_iterator slot=pmtracetot.begin();
	  slot!=pmtracetot.end(); slot++)
	{
	  npesum+=(*slot).second;
	  profpmsum->Fill((*slot).first,(*slot).second);
	  if((*slot).second > npemax) npemax=(*slot).second;
	}
      hpmsumch->Fill(npesum);
      hpmsumpk->Fill(npemax);
       if(npemax!=0)hpmsumatopk->Fill(npesum/npemax);
       
    }
  fCanvasList->Add(new TCanvas("charges","Charges",1000,1000));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2);
  fCanvas->cd(1);
  hpesumch->SetLineColor(kRed);
  hpesumch->Draw();
  fCanvas->cd(2);
  hpmsumch->SetLineColor(kRed);
  hpmsumch->Draw();
  fCanvas->cd(3);
  hadcsumch->SetLineColor(kRed);
  hadcsumch->Draw();
  
  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      hpech[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpech[ipm]->Draw();
      else hpech[ipm]->Draw("same");
      
    }

  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hpmch[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpmch[ipm]->Draw();
      else hpmch[ipm]->Draw("same");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hadcch[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hadcch[ipm]->Draw();
      else hadcch[ipm]->Draw("same");
      
    }
  fCanvasList->Add(new TCanvas("Profiles","Profiles",1000,1000));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2); 

  fCanvas->cd(1);
  profpesum->SetMarkerColor(kRed);
  profpesum->Draw("hist");
  fCanvas->cd(2);
  profpmsum->SetMarkerColor(kRed);
  profpmsum->Draw("hist");
  fCanvas->cd(3);
  profadcsum->SetMarkerColor(kRed);
  profadcsum->Draw("hist");
  
  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      profpe[ipm]->SetLineColor(ipm+1);
       if(ipm==0)profpe[ipm]->Draw("hist");
       else profpe[ipm]->Draw("histsame");
      
    }

  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      profpm[ipm]->SetLineColor(ipm+1);
      if(ipm==0)profpm[ipm]->Draw("hist");
      else profpm[ipm]->Draw("histsame");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      profadc[ipm]->SetLineColor(ipm+1);
      if(ipm==0)profadc[ipm]->Draw("hist");
      else profadc[ipm]->Draw("samehist");
      
    }

 fCanvasList->Add(new TCanvas("Peaks","Peaks  ",1000,1000));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2);
  fCanvas->cd(1);
  hpesumpk->SetLineColor(kRed);
  hpesumpk->Draw();
  fCanvas->cd(2);
  hpmsumpk->SetLineColor(kRed);
  hpmsumpk->Draw();
  fCanvas->cd(3);
  hadcsumpk->SetLineColor(kRed);
  hadcsumpk->Draw();
  
  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hpepk[ipm]->SetLineColor(ipm+1);
       if(ipm==0)hpepk[ipm]->Draw();
       else hpepk[ipm]->Draw("same");
      
    }

  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      hpmpk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpmpk[ipm]->Draw();
      else hpmpk[ipm]->Draw("same");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hadcpk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hadcpk[ipm]->Draw();
      else hadcpk[ipm]->Draw("same");
      
    }
  fCanvasList->Add(new TCanvas("Ratio","Ratio ",1000,1000));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2);

  fCanvas->cd(1);
  hpesumatopk->SetLineColor(kRed);
  hpesumatopk->Draw();
  fCanvas->cd(2);
  hpmsumatopk->SetLineColor(kRed);
  hpmsumatopk->Draw();
  fCanvas->cd(3);
  hadcsumatopk->SetLineColor(kRed);
  hadcsumatopk->Draw();
  
  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      hpeatopk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpeatopk[ipm]->Draw();
      else hpeatopk[ipm]->Draw("same");
      
    }
  
  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hpmatopk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpmatopk[ipm]->Draw();
      else hpmatopk[ipm]->Draw("same");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
    
      hadcatopk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hadcatopk[ipm]->Draw();
      else hadcatopk[ipm]->Draw("same");
      
    }

 

}


void Analyze::DrawMuons()
{
  
  //short adcsum[2][MAXNUMBEROFADCBINS];
  int nmu;
  double npesum,npemax;
  HitStation* sta;
  TProfile *profadcsum;
  TH1F*hadcsumch;
  TH1F*hadcsumpk;
  TH1F*hadcsumatopk;
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);

 
 
  hadcsumch =  new TH1F("adcsumcharge","adcsumcharge",600,0.,600);
  hadcsumpk =  new TH1F("adcsumpeak","adcsumpeak",600,0.,600);
  hadcsumatopk =  new TH1F("adcsumatopk","adcsumatopk",100,0.,10);
  profadcsum = new TProfile("profadcsum","profadcsum",500,0,50,0,50);

  
  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;



  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
     
      sta=&(fStList[0]);
  
     // work with adcsignal  
      npesum=0;
      npemax=0;
      //adcsum =sta->fADC;
      
      // new loop on timebins
      for(int k=0; k<100; k++){
	if(sta->fADC[0][k]!=0){
	  profadcsum->Fill(k,sta->fADC[0][k]);
	  npesum+=sta->fADC[0][k] ;
	  if(sta->fADC[0][k]>npemax) npemax=sta->fADC[0][k];
	
	} 
	else
	  profadcsum->Fill(k,0);
     }// end of loop on k bins
      hadcsumch->Fill(npesum);
      hadcsumpk->Fill(npemax);
      if(npemax!=0)hadcsumatopk->Fill(npesum/npemax);
      
     
    }
  fCanvasList->Add(new TCanvas("muAnalyse"," mu analyse ",800,800));
 
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,2);
  fCanvas->cd(1);
  hadcsumch->SetLineColor(kRed);
  hadcsumch->Draw();
  hadcsumch->Fit("gaus");
  fCanvas->cd(2);
  profadcsum->SetLineColor(kRed);
  profadcsum->Draw();
  fCanvas->cd(3);
  hadcsumpk->SetLineColor(kRed);
  hadcsumpk->Draw();
  hadcsumpk->Fit("gaus");
  fCanvas->cd(4);
  hadcsumatopk->SetLineColor(kRed);
  hadcsumatopk->Draw();
  hadcsumatopk->Fit("gaus");
}

void Analyze::DrawMuons(int opt)
{

  ofstream out;
  out.open("resultfit.out",ios::out);
  // short adcsum[2][MAXNUMBEROFADCBINS];
  map<int,double> trace;
  map<int,double> pmtrace;
  map<int,double> pmtracetot;

  int nmu;
  double npesum,npemax;
  HitStation* sta;

  TH1F*hadcch[NPM];
  TH1F*hadcsumch;
  TH1F*hpech[NPM];
  TH1F*hpesumch;
  TH1F*hpmch[NPM];
  TH1F*hpmsumch;

  TH1F*hadcpk[NPM];
  TH1F*hadcsumpk;
  TH1F*hpepk[NPM];
  TH1F*hpesumpk;
  TH1F*hpmpk[NPM];
  TH1F*hpmsumpk;

  TH1F*hadcatopk[NPM];
  TH1F*hadcsumatopk;
  TH1F*hpeatopk[NPM];
  TH1F*hpesumatopk;
  TH1F*hpmatopk[NPM];
  TH1F*hpmsumatopk;
  
  TProfile *profadcsum;
  TProfile *profpmsum;
  TProfile *profpesum;
  TProfile *profadc[NPM];
  TProfile *profpm[NPM];
  TProfile *profpe[NPM];
  
 
  gStyle->SetOptStat(1001111); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);
 
  hadcsumch =  new TH1F("ADCsumcharge","ADCsumcharge",5000,0.,5000);
  hadcsumpk =  new TH1F("adcsumpeak","adcsumpeak",200,0.,200);
  hadcsumatopk =  new TH1F("adcsumatopk","adcsumatopk",100,0.,10);
  profadcsum = new TProfile("profadcsum","profadcsum",40,0,40,0,450);
  hpesumch =  new TH1F("pesumcharge","pesumcharge",1500,0.,3000);
  hpesumpk =  new TH1F("pesumpeak","pesumpeak",100,0.,100);
  hpesumatopk =  new TH1F("pesumatopk","pesumatopk",100,0.,50);
  profpesum = new TProfile("profpesum","profpesum",1000,0,1000,0,50);
  hpmsumch =  new TH1F("pmsumcharge","pmsumcharge",2000,0.,500);
  hpmsumpk =  new TH1F("pmsumpeak","pmsumpeak",200,0.,1);
  hpmsumatopk =  new TH1F("pmsumatopk","pmsumatopk",200,0.,200);
  profpmsum = new TProfile("profpmsum","profpmsum",1000,0,1000,0,50);

  
   
  
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      char * hname1 = new char[15];
      sprintf(hname1,"ADCcharge_%1d",ipm+1);
      hadcch[ipm] =  new TH1F(hname1,hname1,1000,0.,1000);
      delete[] hname1;
      char* hname2 = new char[15];
      sprintf(hname2,"ADCpeak_%1d",ipm+1);
      hadcpk[ipm] =  new TH1F(hname2,hname2,200,0.,200);
      delete[] hname2;
     char* hname3 = new char[15];
      sprintf(hname3,"ADCatopeak_%1d",ipm+1);
      cout<<hname3<<endl;
      hadcatopk[ipm] =  new TH1F(hname3,hname3,100,0.,10);
       delete[] hname3;
      char* hname4 = new char[15];
      sprintf(hname4,"profADC_%1d",ipm+1);
      profadc[ipm] =  new TProfile(hname4,hname4,40,0.,40,0,50);
      delete[] hname4;
      char* hname5 = new char[15];
      sprintf(hname5,"PEcharge_%1d",ipm+1);
      hpech[ipm] =  new TH1F(hname5,hname5,200,0.,200);
      delete[] hname5;
      char*hname6 = new char[15];
      sprintf(hname6,"PEpeak_%1d",ipm+1);
      hpepk[ipm] =  new TH1F(hname6,hname6,50,0.,50);
      delete[] hname6;
      char*hname7 = new char[15];
      sprintf(hname7,"PEatopeak_%1d",ipm+1);
      hpeatopk[ipm] =  new TH1F(hname7,hname7,100,0.,50);
      delete[] hname7;
      char* hname8 = new char[15];
      sprintf(hname8,"profPE_%1d",ipm+1);
      profpe[ipm] =  new TProfile(hname8,hname8,10000,0.,1000,0,50);
      delete[] hname8;
      char*hname9 = new char[15];
      sprintf(hname9,"pmcharge_%1d",ipm+1);
      hpmch[ipm] =  new TH1F(hname9,hname9,100,0.,20);    
      delete[] hname9;
      char* hname10 = new char[15];
      sprintf(hname10,"pmpeak_%1d",ipm+1);
      hpmpk[ipm] =  new TH1F(hname10,hname10,200,0.,1);
      delete[] hname10;
      char* hname11 = new char[15];
      sprintf(hname11,"pmatopeak_%1d",ipm+1);
      hpmatopk[ipm] =  new TH1F(hname11,hname11,400,0.,200);
      delete[] hname11;
       char* hname12 = new char[15];
      sprintf(hname12,"profpm_%1d",ipm+1);
      profpm[ipm] =  new TProfile(hname12,hname12,1000,0.,1000,0,50);
      delete[] hname12;

    }

  

  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;
  
  
  
  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      sta=&(fStList[0]);
      
      // work with adcsignal  
      npesum=0;
      npemax=0;
      //adcsum =sta->fADC;
       for(int k=0; k<400; k++){
	if(sta->fADC[0][k]!=0){
	  profadcsum->Fill(k,sta->fADC[0][k]);
	  npesum+=sta->fADC[0][k] ;
	  if(sta->fADC[0][k]>npemax) npemax=sta->fADC[0][k];
	  
	} 
	else
	  profadcsum->Fill(k,0);
      }// end of loop on k bins
      if(npesum>0)hadcsumch->Fill(npesum);
      if(npemax>0)hadcsumpk->Fill(npemax);
      if(npemax!=0)hadcsumatopk->Fill(npesum/npemax);
      
      
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  // adcsum =sta->fPMT[ipm].fADC;
	  
	  // new loop on timebins
	  for(int k=0; k<400; k++){
	    if(sta->fPMT[ipm].fADC[0][k]!=0){
	      profadc[ipm]->Fill(k,sta->fPMT[ipm].fADC[0][k]);
	      npesum+=sta->fPMT[ipm].fADC[0][k] ;
	      if(sta->fPMT[ipm].fADC[0][k]>npemax) npemax=sta->fPMT[ipm].fADC[0][k];
	    } 
	    else
	      profadc[ipm]->Fill(k,0);
	  }// end of loop on k bins
	  if(npesum>0)hadcch[ipm]->Fill(npesum);
	  if(npemax>0)hadcpk[ipm]->Fill(npemax);
	  if(npemax!=0)hadcatopk[ipm]->Fill(npesum/npemax);
	  
	}
      
      npesum=0;
      npemax=0;
      
      trace=sta->fTimeProfile;
      
      for(int k=0; k<10000; k++){
	if(trace.count(k)){
	  npesum+=trace[k];
	  profpesum->Fill(k,trace[k]);
	  if(trace[k] > npemax) npemax=trace[k];
	}
	else
	  profpesum->Fill(k,0);
      }
      hpesumch->Fill(npesum);
      hpesumpk->Fill(npemax);
      if(npemax!=0)hpesumatopk->Fill(npesum/npemax);
      
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  
	  trace =sta->fPMT[ipm].fTimeProfile;
	  
	  for(int k=0; k<10000; k++){
	    if(trace.count(k))
	      {
	    
		npesum+=trace[k];
		profpe[ipm]->Fill(k,trace[k]);
		if(trace[k] > npemax) npemax=trace[k];
	      }
	    else
	      profpe[ipm]->Fill(k,0);
	      
	    }
	  hpech[ipm]->Fill(npesum);
	  hpepk[ipm]->Fill(npemax);
	  if(npemax!=0)hpeatopk[ipm]->Fill(npesum/npemax);
	  
	  
	}
      
      
      npesum=0;
      npemax=0;
      
      pmtracetot.clear();
      
      
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  
	  pmtrace =sta->fPMT[ipm].fPMTSignal_hi;
	  
	  for(int k=0; k<10000; k++){
	    if(pmtrace.count(k))
	      {
	      npesum+=pmtrace[k];
	      profpm[ipm]->Fill(k,pmtrace[k]);
	      if(pmtrace[k] > npemax) npemax=pmtrace[k];
	      if( pmtracetot.count(k))
		pmtracetot[k]+=pmtrace[k];
	      else
		pmtracetot[k]+=pmtrace[k];
	      if(pmtrace[k] > npemax) npemax=pmtrace[k];
	      
	    }
	    else
	      {
		profpm[ipm]->Fill(k,0);
	      }
	  }
	  hpmch[ipm]->Fill(npesum);
	  hpmpk[ipm]->Fill(npemax);
	  if(npemax!=0)hpmatopk[ipm]->Fill(npesum/npemax);
	    
	  
	}
      
      npesum=0;
      npemax=0;
      
      for(int k=0; k<10000; k++){
	if(pmtracetot.count(k))
	  {
	    npesum+=pmtracetot[k];
	    profpmsum->Fill(k,pmtracetot[k]);
	    if(pmtracetot[k] > npemax) npemax=pmtracetot[k];
	    
	    
	  }
	else
	  {
	    profpmsum->Fill(k,0);
	  }
      
      }
      for(map<int,double>::const_iterator slot=pmtracetot.begin();
	  slot!=pmtracetot.end(); slot++)
	{
	  npesum+=(*slot).second;
	  profpmsum->Fill((*slot).first,(*slot).second);
	  if((*slot).second > npemax) npemax=(*slot).second;
	}
      hpmsumch->Fill(npesum);
      hpmsumpk->Fill(npemax);
       if(npemax!=0)hpmsumatopk->Fill(npesum/npemax);
       
    }
  fCanvasList->Add(new TCanvas("charges","Charges",1000,1000));
  TF1 *f1 = new TF1("f1","gaus");
  TF1 *f2 = new TF1("f2","expo");
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2);
  fCanvas->cd(1);
  hpesumch->SetLineColor(kRed);
  hpesumch->Draw();
  fCanvas->cd(2);
  hpmsumch->SetLineColor(kRed);
  hpmsumch->Draw();
 fCanvas->cd(3);
  hadcsumch->SetLineColor(kRed);
  hadcsumch->Draw();
   hadcsumch->Fit(f1);
   out<<"ch "<<f1->GetParameter(1)/(int)NPM<<endl;
  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      hpech[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpech[ipm]->Draw();
      else hpech[ipm]->Draw("same");
      
    }

  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hpmch[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpmch[ipm]->Draw();
      else hpmch[ipm]->Draw("same");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hadcch[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hadcch[ipm]->Draw();
      else{
	hadcch[ipm]->Draw("same");
    	hadcch[ipm]->Fit(f1);
	out<<"ch "<<ipm+1<<" "<<f1->GetParameter(1)<<endl;
      }
      //   hadcch[ipm]->Fit("gaus");
    }
  fCanvasList->Add(new TCanvas("Profiles","Profiles",1000,1000));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2); 

  fCanvas->cd(1);
  profpesum->SetMarkerColor(kRed);
  profpesum->Draw("hist");
  fCanvas->cd(2);
  profpmsum->SetMarkerColor(kRed);
  profpmsum->Draw("hist");
  fCanvas->cd(3);
  profadcsum->SetMarkerColor(kRed);
  profadcsum->Draw("hist");
  profadcsum->Fit(f2," "," ",4,25);
  out<<"tau "<<25/f2->GetParameter(1)<<endl;

  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      profpe[ipm]->SetLineColor(ipm+1);
       if(ipm==0)profpe[ipm]->Draw("hist");
       else profpe[ipm]->Draw("histsame");
      
    }

  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      profpm[ipm]->SetLineColor(ipm+1);
      if(ipm==0)profpm[ipm]->Draw("hist");
      else profpm[ipm]->Draw("histsame");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      profadc[ipm]->SetLineColor(ipm+1);
      if(ipm==0)profadc[ipm]->Draw("hist");
      else 
	profadc[ipm]->Draw("samehist");
      
      profadc[ipm]->Fit(f2," "," ",4,25);
   out<<"tau "<<ipm+1<<" "<<25/f2->GetParameter(1)<<endl;

   }
  
 fCanvasList->Add(new TCanvas("Peaks","Peaks  ",1000,1000));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2);
  fCanvas->cd(1);
  hpesumpk->SetLineColor(kRed);
  hpesumpk->Draw();
  fCanvas->cd(2);
  hpmsumpk->SetLineColor(kRed);
  hpmsumpk->Draw();
  fCanvas->cd(3);
  hadcsumpk->SetLineColor(kRed);
  hadcsumpk->Draw();
  hadcsumpk->Fit(f1);
  out<<"pk "<<f1->GetParameter(1)/(int)NPM<<endl;

  
  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hpepk[ipm]->SetLineColor(ipm+1);
       if(ipm==0)hpepk[ipm]->Draw();
       else hpepk[ipm]->Draw("same");
      
    }

  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      hpmpk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpmpk[ipm]->Draw();
      else hpmpk[ipm]->Draw("same");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hadcpk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hadcpk[ipm]->Draw();
      else hadcpk[ipm]->Draw("same");
      hadcpk[ipm]->Fit(f1);
    out<<"pk "<<ipm+1<<" "<<f1->GetParameter(1)<<endl;

    
    }
  fCanvasList->Add(new TCanvas("Ratio","Ratio ",1000,1000));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(3,2);

  fCanvas->cd(1);
  hpesumatopk->SetLineColor(kRed);
  hpesumatopk->Draw();
  fCanvas->cd(2);
  hpmsumatopk->SetLineColor(kRed);
  hpmsumatopk->Draw();
  fCanvas->cd(3);
  hadcsumatopk->SetLineColor(kRed);
  hadcsumatopk->Draw();
  hadcsumatopk->Fit(f1);
  out<<"atop "<<f1->GetParameter(1)<<endl;


  fCanvas->cd(4);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      hpeatopk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpeatopk[ipm]->Draw();
      else hpeatopk[ipm]->Draw("same");
      
    }
  
  fCanvas->cd(5);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
     
      hpmatopk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hpmatopk[ipm]->Draw();
      else hpmatopk[ipm]->Draw("same");
    }
  
  fCanvas->cd(6);
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
    
      hadcatopk[ipm]->SetLineColor(ipm+1);
      if(ipm==0)hadcatopk[ipm]->Draw();
      else hadcatopk[ipm]->Draw("same");
      hadcatopk[ipm]->Fit(f1);
     out<<"atop "<<ipm+1<<" "<<f1->GetParameter(1)<<endl;

 }

  out.close(); 

}
void Analyze::DrawMuons(int n1,int n2,int opt)
{
  
  //short adcsum[2][MAXNUMBEROFADCBINS];
  map<int,double> pmtrace;
 
  HitStation* sta;
  int indix=0;
  int indix2=0;
  TH1F* histo;

  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  int nmu=(int)fTree->GetEntries();

  cout<<nmu<<" muons to draw"<<endl;
  
  if(opt==0){
    /*for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
      for(int j = 0;j<2;j++)
      adcsum[j][i]=0;*/
    
    for(int ne=n1; ne<=n2; ne++)
      {
	fTree->GetEvent(ne);
	if(ne%1000==0) cout<< ne << " particles "<<endl;
	fStList = fEvent->fHitStationList;
	if(fStList.size()!=1) cout<<" Error nb of tanks "<<fStList.size()<<endl;
	sta=&(fStList[0]);
	
	
	if(sta->fNpe>0) 
	  {
	    
	    if ( indix%8==0)
	      {
		
	      indix=0;
	      indix2++;
	      char* canvasname = new char[2];
	      sprintf(canvasname,"c%1d",indix2);
	      
	      fCanvasList->Add(new TCanvas(canvasname,"Traces ",600,800));
	      delete canvasname;
	      fCanvas=(TCanvas *) fCanvasList->Last();
	      fCanvas->Divide(2,4);
	      cout<<"new canvas"<<endl;
	      
	      } 
	    
	    indix++;
	    
	    
	    
	    
	    char*hname = new char[15];
	    
	  
	    sprintf(hname,"SUMADC_%04d",ne);
	    cout<<hname<<endl;
	    histo=  new TH1F(hname,hname,500,0.,500);
	    
	    delete hname;
	    
	    fHisto1->Add(histo);
	    // adcsum=sta->fADC;
	    
	    for(int k=0; k<1000; k++)
	      if(sta->fADC[0][k]!=0){
		histo->Fill(k,sta->fADC[0][k]);
	      }
	    
	    
	    fCanvas->cd(indix);
	    histo->SetLineColor(kRed);
	    histo->Draw();
	    
	    
	  }
	
	
	if( indix%8==0)
	  {
	    fCanvas->Print();
	  }      
	
	
    
      }// end of loop on events
  }//end of opt =0
  if(opt==1){
   
    
    for(int ne=n1; ne<=n2; ne++)
      {
	fTree->GetEvent(ne);
	if(ne%1000==0) cout<< ne << " particles "<<endl;
	fStList = fEvent->fHitStationList;
	if(fStList.size()!=1) cout<<" Error nb of tanks "<<fStList.size()<<endl;
	sta=&(fStList[0]);
	
	
	if(sta->fNpe>0) 
	  {
	    
	    if ( indix%8==0)
	      {
		
	      indix=0;
	      indix2++;
	      char* canvasname = new char[2];
	      sprintf(canvasname,"c%1d",indix2);
	      
	      fCanvasList->Add(new TCanvas(canvasname,"Traces ",600,800));
	      delete canvasname;
	      fCanvas=(TCanvas *) fCanvasList->Last();
	      fCanvas->Divide(2,4);
	      cout<<"new canvas"<<endl;
	      
	      } 
	    
	    indix++;
	    
	    
	    
	    
	    char*hname = new char[15];
	    
	  
	    sprintf(hname,"SUMPMT_%04d",ne);
	    cout<<hname<<endl;
	    histo=  new TH1F(hname,hname,500,0.,500);
	    
	    delete hname;
	    
	    for(int ipm=0;ipm<NPM;ipm++){
	      pmtrace=sta->fPMT[ipm].fPMTSignal_hi;   
	      
	      for(map<int,double>::const_iterator slot=pmtrace.begin();
		  slot!=pmtrace.end(); slot++)
		
		histo->Fill((*slot).first,(*slot).second);
		
	    }
	    
	    fCanvas->cd(indix);
	    histo->SetLineColor(kRed);
	    histo->Draw();
	    
	    
	  }
	
	
	if( indix%8==0)
	  {
	    fCanvas->Print();
	  }      
	
	
	
      }// end of loop on events
  }//end of opt =1
  
}

void Analyze::DrawMuonsCharge(int ipm)
{
  
  //short adcsum[2][MAXNUMBEROFADCBINS];
  

  int nmu;
  double npesum=0;
  double npemax=0;
  HitStation* sta;
 
  TH1F*hadcch;
  TH1F*hadcsumch;
  TH1F*hadcpk;
  TH1F*hadcsumpk;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLineWidth(2);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);

  char*hname ;
  
  hadcsumch =  new TH1F("ADCsumcharge","ADCsumcharge",1500,0.,1500);
  hname = new char[11];
  sprintf(hname,"ADCcharge_pm%1d",ipm);
  hadcch =  new TH1F(hname,hname,1500,0.,1500);
  delete hname;
   
  hadcsumpk =  new TH1F("ADCsumpeak","ADCsumpeak",500,0.,500);
  hname = new char[11];
  sprintf(hname,"ADCpeak_pm%1d",ipm);
  hadcpk =  new TH1F(hname,hname,500,0.,500);
  delete hname;
  
  
  
  
  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
     for(int j = 0;j<2;j++)
     adcsum[j][i]=0;*/
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;

  
  
  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      sta=&(fStList[0]);
  
      // work with adcsignal  
      npesum=0;
      npemax=0;
      
      //adcsum =sta->fADC;
      
      for(int k=0; k<100; k++)
	if(sta->fADC[0][k]!=0){
	  npesum+=sta->fADC[0][k] ;
	  if(sta->fADC[0][k]>npemax)npemax=sta->fADC[0][k] ;
	}
      
      hadcsumch->Fill(npesum/(float)NPM);
      hadcsumpk->Fill(npemax/(float)NPM);
      
      
      npesum=0;
      npemax=0;
      
      //adcsum =sta->fPMT[ipm-1].fADC;
      
      // new loop on timebins
      for(int k=0; k<100; k++)
	if(sta->fPMT[ipm-1].fADC[0][k]!=0){
	  
	  npesum+=sta->fPMT[ipm-1].fADC[0][k] ;
	  if(sta->fPMT[ipm-1].fADC[0][k]>npemax)npemax=sta->fPMT[ipm-1].fADC[0][k] ;
	}
      hadcch->Fill(npesum);
      hadcpk->Fill(npemax);
      
      
      
     
      
    }
  
  fCanvasList->Add(new TCanvas("ADC charge"," ADC Charge ",800,400));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,1);
  fCanvas->cd(2);
  hadcsumch->SetLineColor(1);
  hadcsumch->Draw();
  hadcch->SetLineColor(kRed);
  hadcch->Draw("same");
  fCanvas->cd(1);
  hadcsumpk->SetLineColor(1);
  hadcsumpk->Draw();
  hadcpk->SetLineColor(kRed);
  hadcpk->Draw("same");
 


}


void Analyze::DrawMuonsFE()
{
  
   HitStation* sta;

  int nmu;
  double sumFE;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);

  //  double fesum[MAXNUMBEROFTIMEBINS];
  
  TH1F* hcharge;

  int nmax=MAXNUMBEROFTIMEBINS ;
  
  //  for(int i = 0;i<nmax;i++)
  // fesum[i]=0;
  
  
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;

  hcharge= new TH1F("chargeFE","chargeFE",100,0,50);
  fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
 

  for(int ne=0; ne<nmu; ne++)
    {
      
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
 
      // for(int i = 0;i<nmax;i++) fesum[i]=0;
 
      sta=&(fStList[0]);
      sumFE=0;

      for (int ipm=0;ipm<NPM;ipm++)
	{
	  // fesum=sta->fPMT[ipm].fFEFilter_hi;
    
	  for(int k=0; k<nmax; k++)
	    if(sta->fPMT[ipm].fFEFilter_hi[k]!=0){
	      sumFE+=sta->fPMT[ipm].fFEFilter_hi[k];
	    }
	}
      cout<<sumFE<<endl;
      hcharge->Fill(sumFE);

    }

  hcharge->Draw();
}

void Analyze::DrawMuonSpectrum()
{
  FILE *fp = fopen("mudecay.out","w");
  FILE *fp1 = fopen("mudecay_1.out","w");
  FILE *fp2 = fopen("mudecay_2.out","w");
  FILE *fp3 = fopen("mudecay_3.out","w");
 
  ofstream out;
  out.open("decaytest.out",ios::out);
  //short adcsum[2][MAXNUMBEROFADCBINS];
  //short adc1[2][MAXNUMBEROFADCBINS];
  //short adc2[2][MAXNUMBEROFADCBINS];
  //short adc3[2][MAXNUMBEROFADCBINS];
  int nmu;
  double mutime0;
  double mutime;
  double mutime1;
  double dmutime;
  int th=3;
  HitStation* sta;
  
  TH1F* htime;
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);
  
  //out.open("mudecay.out",ios::out);
  
  htime =  new TH1F("time","time",10000,0.,10000);
 
  
  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;
  
  int intpk0=0; 
  
  for(int ne=0; ne<nmu; ne++)
    {
      int binpk1=0;
      int binpk2=0;
      int endpk1=0;
      int endpk2=0;
      int debpk1=0;
      int debpk2=0;
      int intpk1=0;
      int intpk2=0;
      int intpk1_1=0;
      int intpk2_1=0;
      int intpk1_2=0;
      int intpk2_2=0;
      int intpk1_3=0;
      int intpk2_3=0;
   
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      sta=&(fStList[0]);
      // adcsum =sta->fADC;
      //adc1 =sta->fPMT[0].fADC;
      //adc2 =sta->fPMT[1].fADC;
      //adc3 =sta->fPMT[2].fADC;
      
     if(ne>0){
	
	mutime= sta->fTime0/25.;
	for(int k=0; k<400; k++){
	  if(sta->fPMT[0].fADC[0][k]>th && sta->fPMT[1].fADC[0][k]>th && sta->fPMT[2].fADC[0][k]>th && debpk1==0){
	    debpk1=1;
	    binpk1=k;
	  }
	  if(debpk1==1 && sta->fADC[0][k]<1)endpk1++;
	  if(debpk1==1 && endpk1<3){
	    intpk1+=sta->fADC[0][k] ;
	    intpk1_1+=sta->fPMT[0].fADC[0][k] ;
	    intpk1_2+=sta->fPMT[1].fADC[0][k] ;
	    intpk1_3+=sta->fPMT[2].fADC[0][k] ;
	  }
	  if(endpk1>3){
	    if(sta->fPMT[0].fADC[0][k]>th && sta->fPMT[1].fADC[0][k]>th && sta->fPMT[2].fADC[0][k]>th && debpk2==0){
	      debpk2=1;
	      binpk2=k;
	    }
	    if(debpk2==1 && sta->fADC[0][k]<1)endpk2++;
	    if(debpk2==1 && endpk2<3){
	      intpk2+=sta->fADC[0][k] ;
	      intpk2_1+=sta->fPMT[0].fADC[0][k] ;
	      intpk2_2+=sta->fPMT[1].fADC[0][k] ;
	      intpk2_3+=sta->fPMT[2].fADC[0][k] ;
	    }
	    if(endpk2>3)break;
	  }
	}
	if(intpk1>0){
	mutime1=mutime+binpk1;
	dmutime=mutime1-mutime0;
	htime->Fill(dmutime);
	fprintf(fp,"%.10e %d\n",mutime1,intpk1);
	fprintf(fp1,"%.10e %d\n",mutime1,intpk1_1);
	fprintf(fp2,"%.10e %d\n",mutime1,intpk1_2);
	fprintf(fp3,"%.10e %d\n",mutime1,intpk1_3);
//out<<mutime1<<" "<<intpk1<<endl;
	mutime0=mutime1;
	intpk0=intpk1;
	}
	

	//if two peaks in the trace
	
	if(intpk2>0){
	  mutime1=mutime+binpk2;
	  dmutime=mutime1-mutime0;
	  htime->Fill(dmutime);
	  //out<<ne<<" "<<mutime0<<" "<<mutime1<<" "<<intpk1<<" "<<intpk2<<endl;
	  mutime0=mutime1;
	  //  out<<mutime1<<" "<<intpk2<<endl;
	  fprintf(fp,"%.10e %d\n",mutime1,intpk2);
	  fprintf(fp1,"%.10e %d\n",mutime1,intpk2_1);
	  fprintf(fp2,"%.10e %d\n",mutime1,intpk2_2);
	  fprintf(fp3,"%.10e %d\n",mutime1,intpk2_3);
	  intpk0=intpk2;
	}
	
	out<<ne<<" "<<binpk1<<" "<<binpk2<<" "<<intpk1<<" "<<intpk2<<" "<<intpk1_1<<" "<<intpk2_1<<" "<<intpk1_2<<" "<<intpk2_2<<" "<<intpk1_3<<" "<<intpk2_3<<endl;
     }
     else
       {
	 mutime0=sta->fTime0;
	 for(int k=0; k<400; k++){
	   if(sta->fADC[0][k]>0){
	     intpk0+=sta->fADC[0][k];
	   }
	    
	  }
	}
     // dmutime= sta->fTime0;
     // work with adcsignal  
     // cout<<dmutime<<endl;
     
      
      
      
      
      
    
      
    }
  fCanvasList->Add(new TCanvas("mutime"," mu time ",600,600));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  htime->SetLineColor(kRed);
  htime->Draw();
}


void Analyze::DrawMuonTimes()
{
  ifstream files;
  char * filename;
  TFile * file;
  filename = new char[40];
  int nfile,nevent,nsta,nmu,id,mumax;
  double dist,tmu,tmusc;
  HitStation* sta;
  TH1F *hmutime[2];
  TH1F *hmutimesc[2];
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);

  for (int i=0;i<2;i++){
    hmutime[i]= new TH1F("mutime","mutime",50,0,5000);
    hmutimesc[i]= new TH1F("mutimesc","mutimesc",50,0,5000);
  }
 
  files.open("files.inp",ios::in);
  files >> mumax;
  files>>nfile;
  for (int i=0;i<nfile;i++)
    {
      files>>filename;
      
      file = 0;
      fTree=0;
      fBranch=0;
      fEvent=0;
      
      
      file = new TFile(filename,"READ");	 // Open file	
      file->ls();
     
      fTree = (TTree *) file->Get("TreeEvent");      // Get Tree
      
      fTree->Print(); 
      
      

      fBranch= fTree->GetBranch("event");
      fTree->SetBranchAddress("event",&fEvent);
      
      cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
      
    
      nevent=(int)fTree->GetEntries();
      cout<<nevent<<" showers  to analyze"<<endl;

 
      
      for(int ne=0; ne<nevent; ne++)
	{
	  fTree->GetEvent(ne);
	  if(ne%1000==0) cout<< ne << " particles "<<endl;
	  fStList = fEvent->fHitStationList;
	  id =  fEvent->fPrimary;
	  nsta = fStList.size();
	  cout<<id<<" "<<nsta<<endl;
	  for (int ns=0;ns<nsta;ns++)
	    {
	      sta=&(fStList[ns]);
	      nmu=sta->fNmu;
	      dist = sta->fR_sf;
	      cout<<nmu<<"\t  muons at  "<<dist<<" m"<<endl;
	      if (nmu<mumax)
		for (int nm=0;nm<nmu;nm++){
		  tmu = sta->fMuTimes[nm];
		  tmusc = tmu *1000 /dist;	
		  if(id==31 || id ==-31)
		    {
		      hmutime[0]->Fill(tmu,1);
		      hmutimesc[0]->Fill(tmusc,1);
		    }
		  if(id==2656 )
		    {
		      hmutime[1]->Fill(tmu,1);
		      hmutimesc[1]->Fill(tmusc,1);
		    }
		}
	    }
	  
	}//end of loop en event number
      file->Close(0);
    }//end of loop on files
  files.close();
  
  fCanvasList->Add(new TCanvas("mutimes"," mutimes ",600,600));
  
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(1,2);
  fCanvas->cd(1);
  
  hmutime[0]->SetMarkerStyle(20);
  hmutime[0]->SetMarkerColor(kRed);
  hmutime[0]->Draw();
  hmutime[1]->SetMarkerStyle(21);
  hmutime[1]->SetMarkerColor(kBlue);
  hmutime[1]->Draw("same");
  fCanvas->cd(2);
  hmutimesc[0]->SetMarkerStyle(20);
  hmutimesc[0]->SetMarkerColor(kRed);
  hmutimesc[0]->Draw();
  hmutimesc[1]->SetMarkerStyle(21);
  hmutimesc[1]->SetMarkerColor(kBlue);
  hmutimesc[1]->Draw("same");
}

void Analyze::DrawOneFE(int num,int ipm)
{
  // double fesum[MAXNUMBEROFTIMEBINS];
  
  TH1F* histo;
  
  int nmax=MAXNUMBEROFTIMEBINS ;
  
  // for(int i = 0;i<nmax;i++)
  // fesum[i]=0;
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  
  int nmu=fStList[num].fNmu;
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
 
  ipm--;
  
     
      fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      
      sprintf(hname,"FE%1d_%04d_%04d",ipm+1,id,idistance);
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,5000,0,5000);
      delete hname;
      
      fHisto1->Add(histo);
     


      //      fesum=fStList[num].fPMT[ipm].fFEFilter_hi;
      for(int k=0; k<nmax; k++)
	if(fStList[num].fPMT[ipm].fFEFilter_hi[k]!=0){
	  histo->Fill(k,fStList[num].fPMT[ipm].fFEFilter_hi[k]);
	}
    
     
      histo->SetLineColor(kRed);
      histo->Draw();
      
	
   
}


void Analyze::DrawPulseShape()
{
  map<int,double> pmtrace;
  TProfile *profpm[NPM];
  int nmu;
  double npesum,npemax;
  HitStation* sta;

  gStyle->SetOptStat(1001111); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);

 
  
   
  
  
  for (int ipm=0;ipm<NPM;ipm++)
    {
      
      char* hname12 = new char[15];
      sprintf(hname12,"profpm_%1d",ipm+1);
      profpm[ipm] =  new TProfile(hname12,hname12,1000,0.,1000,0,50);
      delete[] hname12;

    }

  

  
  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;
  
  
  
  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      sta=&(fStList[0]);
      
      // work with adcsignal  
   
      npesum=0;
      npemax=0;
      
   
      for (int ipm=0;ipm<NPM;ipm++)
	{
	  npesum=0;
	  npemax=0;
	  
	  pmtrace =sta->fPMT[ipm].fPMTSignal_hi;
	 
	  for(map<int,double>::const_iterator slot=pmtrace.begin();
	      slot!=pmtrace.end(); slot++)
	    {
	      npesum+=(*slot).second;
	      profpm[ipm]->Fill(64+(*slot).first,(*slot).second);
	     
	      
	    }
	  
	  
	}
    }
  

  fCanvasList->Add(new TCanvas("Profiles","Profiles",800,600));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
 
  profpm[0]->SetLineColor(4);
  profpm[0]->Draw();
  
}

void Analyze::DrawOneFE(int num,int ipm,int ichan)
{
  //double fesum[MAXNUMBEROFTIMEBINS];
  
  TH1F* histo;
  
  int nmax=MAXNUMBEROFTIMEBINS ;
  
  // for(int i = 0;i<nmax;i++)
  //  fesum[i]=0;
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  
  int nmu=fStList[num].fNmu;
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
 
  ipm--;
  
     
      fCanvasList->Add(new TCanvas("c1","Trace ",600,800));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      
      sprintf(hname,"FE%1d_%04d_%04d",ipm+1,id,idistance);
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,5000,0,5000);
      delete hname;
      
      fHisto1->Add(histo);
     


      // if(ichan==0)fesum=fStList[num].fPMT[ipm].fFEFilter_hi;
      //else fesum=fStList[num].fPMT[ipm].fFEFilter_lo;
      if(ichan==0)
	for(int k=0; k<nmax; k++)
	  if(fStList[num].fPMT[ipm].fFEFilter_hi[k]!=0){
	    histo->Fill(k,fStList[num].fPMT[ipm].fFEFilter_hi[k]);
	  }
      
	  else 
	    for(int k=0; k<nmax; k++)
	      if(fStList[num].fPMT[ipm].fFEFilter_lo[k]!=0){
		histo->Fill(k,fStList[num].fPMT[ipm].fFEFilter_lo[k]);
	  }  
 
      histo->SetLineColor(kRed);
      histo->Draw();
      
	
   
}





void Analyze::DrawOnePM(int num,int ipm)
{
  map<int,double> trace;
  TH1F* histo;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

 
  ipm--;
  
  fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  char*hname = new char[15];
  int idistance = (int)fStList[num].fR_sf;
  int id = fStList[num].fId;
  
  sprintf(hname,"PMT%1d_%04d_%04d",ipm+1,id,idistance);
  
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,5000,0.,5000);
  delete hname;
      
  fHisto1->Add(histo);
  trace=fStList[num].fPMT[ipm].fPMTSignal_hi;
  
  for(map<int,double>::const_iterator slot=trace.begin();
      slot!=trace.end(); slot++)
    {
      if(((*slot).first)<10000)
	histo->Fill(((*slot).first),(*slot).second);
	  
    }
  
  
  
  histo->SetLineColor(kRed);
  histo->Draw();
  
	
}

void Analyze::DrawOnePM(int num,int ipm,int ichan)
{
  map<int,double> trace;
  TH1F* histo;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

 
  ipm--;
  
  fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  char*hname = new char[15];
  int idistance = (int)fStList[num].fR_sf;
  int id = fStList[num].fId;
  
  sprintf(hname,"PMT%1d_%04d_%04d",ipm+1,id,idistance);
  
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,5000,0.,5000);
  delete hname;
      
  fHisto1->Add(histo);
  
  if(ichan==0)trace=fStList[num].fPMT[ipm].fPMTSignal_hi;
  else if (ichan==0)trace=fStList[num].fPMT[ipm].fPMTSignal_lo;
  else {
    cout<<" ichan not valide"<<endl;
    return;
  }

  
  for(map<int,double>::const_iterator slot=trace.begin();
      slot!=trace.end(); slot++)
    {
      if(((*slot).first)<10000)
	histo->Fill(((*slot).first),(*slot).second);
	  
    }
  
  
  
  histo->SetLineColor(kRed);
  histo->Draw();
  
	
}

void Analyze::DrawOneProfile(int num)
{
  map<int,double> trace;
  map<int,double> trace_mu;
  TH1F* histo;
  TH1F* histo_mu;
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  int nmu=fStList[num].fNmu;
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
  if(fStList[num].fNpe>0) 
    {
      fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      sprintf(hname,"SUMPMT_%04d_%04d",id,idistance);
      
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,5000,0.,5000);
      histo_mu=  new TH1F(hname,hname,5000,0.,5000);
      delete hname;
      
      fHisto1->Add(histo);
      fHisto1->Add(histo_mu);
      trace=fStList[num].fTimeProfile;
      
      for(map<int,double>::const_iterator slot=trace.begin();
	  slot!=trace.end(); slot++)
	{
	  if(((*slot).first)<10000)
	    histo->Fill(((*slot).first),(*slot).second);
	  
	}
      
      
      trace_mu=fStList[num].fTimeProfile_mu;
      
      for(map<int,double>::const_iterator slot=trace_mu.begin();
	  slot!=trace_mu.end(); slot++)
	{
	  if(((*slot).first)<10000)
	    histo_mu->Fill(((*slot).first),(*slot).second);
	}
      
      histo->SetLineColor(kRed);
      histo_mu->SetLineColor(kCyan);
      histo->Draw();
      histo_mu->Draw("SAME");
      
      
    }
  
}

void Analyze::DrawOneProfile(int num,int ipm)
{
  map<int,double> trace;
  map<int,double> trace_mu;
  TH1F* histo;
  TH1F* histo_mu;
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  int nmu=fStList[num].fNmu;

  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }

  ipm--;
  
  if(fStList[num].fNpe>0) 
    {
      fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      
      sprintf(hname,"ADC%1d_%04d_%04d",ipm+1,id,idistance);
      
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,5000,0.,5000);
      histo_mu=  new TH1F(hname,hname,5000,0.,5000);
      delete hname;
      
      fHisto1->Add(histo);
      fHisto1->Add(histo_mu);
      trace=fStList[num].fPMT[ipm].fTimeProfile;
	  
      for(map<int,double>::const_iterator slot=trace.begin();
	  slot!=trace.end(); slot++)
	{
	  if(((*slot).first)<10000)
	    histo->Fill(((*slot).first),(*slot).second);
	  
	}

      
      trace_mu=fStList[num].fPMT[ipm].fTimeProfile_mu;
      
      for(map<int,double>::const_iterator slot=trace_mu.begin();
	  slot!=trace_mu.end(); slot++)
	{
	  if(((*slot).first)<10000)
	    histo_mu->Fill(((*slot).first),(*slot).second);
	}
      
      histo->SetLineColor(kRed);
      histo_mu->SetLineColor(kCyan);
      histo->Draw();
      histo_mu->Draw("SAME");
      
	      
    }
     
}


void Analyze::DrawOneTrace(int num)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  TH1F* histo;
  TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS;

  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
     for(int j = 0;j<2;j++)
     adcsum[j][i]=0;*/

  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  int nmu=fStList[num].fMuTimes.size();
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
 
  fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  char*hname = new char[15];
  int idistance = (int)fStList[num].fR_sf;
  int id = fStList[num].fId;
  sprintf(hname,"SUMADC_%04d_%04d",id,idistance);
  
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,400,0.,400);
  histo_mu=  new TH1F(hname,hname,400,0.,400);
  delete hname;
  
  fHisto1->Add(histo);
  fHisto1->Add(histo_mu);
  
  //adcsum=fStList[num].fADC;
  for(int k=0; k<nmax; k++)
    if(fStList[num].fADC[0][k]!=0){
      histo->Fill(k,fStList[num].fADC[0][k]);
    }
  
  //adcsum=fStList[num].fADC_mu;
  for(int k=0; k<nmax; k++)
    if(fStList[num].fADC_mu[0][k]!=0){
      histo_mu->Fill(k,fStList[num].fADC_mu[0][k]);
    }
  
  histo->SetLineColor(kRed);
  histo_mu->SetLineColor(kCyan);
  histo->Draw();
  histo_mu->Draw("SAME");
  
  
  
}

void Analyze::DrawOneTrace(int num,int ipm)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  
  TH1F* histo;
  TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS ;
  
  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  
  int nmu=fStList[num].fNmu;
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
 
  ipm--;
  
     
      fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      
      sprintf(hname,"ADC%1d_%04d_%04d",ipm+1,id,idistance);
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,400,0,400);
      histo_mu=  new TH1F(hname,hname,400,0,400);
      delete hname;
      
      fHisto1->Add(histo);
      fHisto1->Add(histo_mu);



      //      adcsum=fStList[num].fPMT[ipm].fADC;
      for(int k=0; k<nmax; k++)
	if(fStList[num].fPMT[ipm].fADC[0][k]!=0){
	  histo->Fill(k,fStList[num].fPMT[ipm].fADC[0][k]);
	}
    
      // adcsum=fStList[num].fPMT[ipm].fADC_mu;
      for(int k=0; k<nmax; k++)
	if(fStList[num].fPMT[ipm].fADC_mu[0][k]!=0){
	  histo_mu->Fill(k,fStList[num].fPMT[ipm].fADC_mu[0][k]);
	}
     
      histo->SetLineColor(kRed);
      histo_mu->SetLineColor(kCyan);
      histo->Draw();
      histo_mu->Draw("SAME");
      
	
   
}
void Analyze::DrawOneTraceInVEM(int num,int ipm,int ichan)
{
  //short adcsum[2][MAXNUMBEROFADCBINS];
  
  TH1F* histo;
  // TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS ;
  
  // for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
  //for(int j = 0;j<2;j++)
  //  adcsum[j][i]=0;
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  
  int nmu=fStList[num].fNmu;
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
 
  ipm--;
  
     
      fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      
      sprintf(hname,"ADC%1d_%04d_%04d",ipm+1,id,idistance);
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,400,0,400);
      //     histo_mu=  new TH1F(hname,hname,400,0,400);
      delete hname;
      
      fHisto1->Add(histo);
 //      fHisto1->Add(histo_mu);



      // adcsum=fStList[num].fPMT[ipm].fADC;
      cout<<VEMPEAKVALUEINADC<<endl;
      for(int k=0; k<nmax; k++)
	if(fStList[num].fPMT[ipm].fADC[ichan][k]!=0){
	  histo->Fill(k,(fStList[num].fPMT[ipm].fADC[ichan][k])/VEMPEAKVALUEINADC*32.);
	}
    
      // adcsum=fStList[num].fPMT[ipm].fADC_mu;
//       for(int k=0; k<nmax; k++)
// 	if(fStList[num].fPMT[ipm].fADC_mu[ichan][k]!=0){
// 	  histo_mu->Fill(k,fStList[num].fPMT[ipm].fADC_mu[ichan][k]);
// 	}
     
      histo->SetLineColor(kRed);
      //      histo_mu->SetLineColor(kCyan);
      histo->Draw();
//       histo_mu->Draw("SAME");
      
	
   
}
void Analyze::DrawOneTrace(int num,int ipm,int ichan)
{
  //short adcsum[2][MAXNUMBEROFADCBINS];
  
  TH1F* histo;
  TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS ;
  
  // for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
  //for(int j = 0;j<2;j++)
  //  adcsum[j][i]=0;
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  
  int nmu=fStList[num].fNmu;
  if(nmu<20) 
    {
      for (int nm=0;nm<nmu;nm++)
	{
	  cout<<"muon "<< nm<<"\t time = "<<fStList[num].fMuTimes[nm]<<endl;
	}
    }
 
  ipm--;
  
     
      fCanvasList->Add(new TCanvas("c1","Trace ",400,300));
      fCanvas=(TCanvas *) fCanvasList->Last();
      
      char*hname = new char[15];
      int idistance = (int)fStList[num].fR_sf;
      int id = fStList[num].fId;
      
      sprintf(hname,"ADC%1d_%04d_%04d",ipm+1,id,idistance);
      cout<<hname<<endl;
      histo=  new TH1F(hname,hname,400,0,400);
      //     histo_mu=  new TH1F(hname,hname,400,0,400);
      delete hname;
      
      fHisto1->Add(histo);
 //      fHisto1->Add(histo_mu);



      // adcsum=fStList[num].fPMT[ipm].fADC;
      for(int k=0; k<nmax; k++)
	if(fStList[num].fPMT[ipm].fADC[ichan][k]!=0){
	  histo->Fill(k,fStList[num].fPMT[ipm].fADC[ichan][k]);
	}
    
      // adcsum=fStList[num].fPMT[ipm].fADC_mu;
//       for(int k=0; k<nmax; k++)
// 	if(fStList[num].fPMT[ipm].fADC_mu[ichan][k]!=0){
// 	  histo_mu->Fill(k,fStList[num].fPMT[ipm].fADC_mu[ichan][k]);
// 	}
     
      histo->SetLineColor(kRed);
      //      histo_mu->SetLineColor(kCyan);
      histo->Draw();
//       histo_mu->Draw("SAME");
      
	
   
}

void Analyze::DrawTraces()
{
  TH1F* histo;
  // TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS;
  
  //short adcsum[2][MAXNUMBEROFADCBINS];
  
  /*  for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  int indix=0;
  int indix2=0;
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  vector<HitStation>::iterator sta;

  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
    {
      // if(sta->fNpe>10) {
        if(sta->fT2ToT||sta->fT1Threshold) 
  	{
	  
	      if ( indix%8==0)
		{
		  
		  indix=0;
		  indix2++;
		  char* canvasname = new char[2];
		  sprintf(canvasname,"c%1d",indix2);
		 
		  fCanvasList->Add(new TCanvas(canvasname,"Traces ",600,800));
		  delete canvasname;
		  fCanvas=(TCanvas *) fCanvasList->Last();
		  fCanvas->Divide(2,4);
		  cout<<"new canvas"<<endl;
		  
		} 
	   
	      indix++;
	      
	 
	      
	      
	      char*hname = new char[15];
	      int idistance = (int)sta->fR_sf;
	      int id = sta->fId;

	      
	      sprintf(hname,"SUMADC_%04d_%04d",id,idistance);
	      cout<<hname<<endl;
	      histo=  new TH1F(hname,hname,400,0.,400);
	      //  histo_mu=  new TH1F(hname,hname,400,0.,400);
	      delete hname;
	      
	      fHisto1->Add(histo);
	      // fHisto1->Add(histo_mu);
	      
	      // adcsum=sta->fADC;
	      for(int k=0; k<nmax; k++)
		if(sta->fADC[0][k]!=0){
		  histo->Fill(k,sta->fADC[0][k]);
		  cout<<id<<" "<<k<<" "<<sta->fADC[0][k]<<endl;
		}
	  
	      //adcsum=sta->fADC_mu;
	     //  for(int k=0; k<nmax; k++)
// 		if(sta->fADC_mu[0][k]!=0){
// 		  histo_mu->Fill(k,sta->fADC_mu[0][k]);
// 		}
	      
	 
	      fCanvas->cd(indix);
	      histo->SetLineColor(kRed);
	      //  histo_mu->SetLineColor(kCyan);
	      histo->Draw();
	      //    histo_mu->Draw("SAME");
	      
	
      
      
	      // if( indix%8==0)
// 		{
// 		  fCanvas->Print();
// 		}      
	      
	
	}// end of loop on station
  }
}

void Analyze::DrawTraces(int ipm)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  TH1F* histo;
  TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS;

  /*  for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/

 
  int indix=0;
  int indix2=0;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  vector<HitStation>::iterator sta;

  ipm--;
  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
    {
     
      if(sta->fNpe>0) 
	{
	  
	      if ( indix%8==0)
		{
		  
		  indix=0;
		  indix2++;
		  char* canvasname = new char[2];
		  sprintf(canvasname,"c%1d",indix2);
		  
		  fCanvasList->Add(new TCanvas(canvasname,"Traces ",600,800));
		  delete canvasname;
		  fCanvas=(TCanvas *) fCanvasList->Last();
		  fCanvas->Divide(2,4);
		  cout<<"new canvas"<<endl;
		  
		} 
	      
	      indix++;
	      
	      
	      
	      
	      char*hname = new char[15];
	      int idistance = (int)sta->fR_sf;
	      
	      int id = sta->fId;
	      
	      sprintf(hname,"PMT%1d_%04d_%04d",ipm+1,id,idistance);
	      
	      
	      cout<<hname<<endl;
	      histo=  new TH1F(hname,hname,400,0.,400);
	      histo_mu=  new TH1F(hname,hname,400,0.,400);
	      delete hname;
	      
	      fHisto1->Add(histo);
	      fHisto1->Add(histo_mu);
	      
	      
	      // adcsum=sta->fPMT[ipm].fADC;
	      for(int k=0; k<nmax; k++)
		if(sta->fPMT[ipm].fADC[0][k]!=0){
		  histo->Fill(k,sta->fPMT[ipm].fADC[0][k]);
		}
	      
	      //adcsum=sta->fPMT[ipm].fADC_mu;
	      for(int k=0; k<nmax; k++)
		if(sta->fPMT[ipm].fADC_mu[0][k]!=0){
		  histo_mu->Fill(k,sta->fPMT[ipm].fADC_mu[0][k]);
		}
	      fCanvas->cd(indix);
	      histo->SetLineColor(kRed);
	      histo_mu->SetLineColor(kCyan);
	      histo->Draw();
	      histo_mu->Draw("SAME");
	      
	      
	}
      
      
      if( indix%8==0)
	{
	  fCanvas->Print();
	}      
      
    }// end of loop on station
  
}



void Analyze::DrawTracesInVem()
{
  TH1F* histo;
  TH1F* histo_mu;
  
  int nmax=MAXNUMBEROFADCBINS;
  
  // short adcsum[2][MAXNUMBEROFADCBINS];
  
  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/
  
  int indix=0;
  int indix2=0;
  
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  vector<HitStation>::iterator sta;

  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
    {
      double sum=0;
      if(sta->fNpe>0 && (sta->fT1Threshold || sta->fT2ToT)) 
	{
	 
	      if ( indix%8==0)
		{
		  
		  indix=0;
		  indix2++;
		  char* canvasname = new char[2];
		  sprintf(canvasname,"c%1d",indix2);
		 
		  fCanvasList->Add(new TCanvas(canvasname,"Traces ",600,800));
		  delete canvasname;
		  fCanvas=(TCanvas *) fCanvasList->Last();
		  fCanvas->Divide(2,4);
		  cout<<"new canvas"<<endl;
		  
		} 
	   
	      indix++;
	      
	 
	      
	      
	      char*hname = new char[15];
	      int idistance = (int)sta->fR_sf;
	      int id = sta->fId;

	      
	      sprintf(hname,"SUMADC_%04d_%04d",id,idistance);
	      cout<<hname<<endl;
	      histo=  new TH1F(hname,hname,nmax,0.,nmax);
	      histo_mu=  new TH1F(hname,hname,nmax,0.,nmax);
	      delete hname;
	      
	      fHisto1->Add(histo);
	      fHisto1->Add(histo_mu);
	      
	      //adcsum=sta->fADC;
	      for(int k=0; k<nmax; k++)
		if(sta->fADC[0][k]!=0){
		  histo->Fill(k,sta->fADC[0][k]/VEMPEAKVALUEINADC);
		  sum+=sta->fADC[0][k];
		}
	  
	      //adcsum=sta->fADC_mu;
	      for(int k=0; k<nmax; k++)
		if(sta->fADC_mu[0][k]!=0){
		  histo_mu->Fill(k,sta->fADC_mu[0][k]/VEMPEAKVALUEINADC);
		}
	      
	 
	      fCanvas->cd(indix);
	      histo->SetLineColor(kRed);
	      histo_mu->SetLineColor(kCyan);
	      histo->Draw();
	      histo_mu->Draw("SAME");
	      
	      cout<<"signal of station "<< sta->fId<<" = "<<sum/VEMCHARGEVALUEINADC/(int)NPM<<" VEM"<<endl;

      
      
	      if( indix%8==0)
		{
		  fCanvas->Print();
		}      
	      
	
	}// end of loop on station
    }
}


void Analyze::DrawTracesInVem(int ipm)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  TH1F* histo;
  TH1F* histo_mu;
  int nmax=MAXNUMBEROFADCBINS;

  /*  for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    for(int j = 0;j<2;j++)
    adcsum[j][i]=0;*/

 
  int indix=0;
  int indix2=0;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  vector<HitStation>::iterator sta;

  ipm--;
  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
    {
      double sum=0;
      if(sta->fNpe>0 &&(sta->fT1Threshold || sta->fT2ToT)) 
	{
	  
	      if ( indix%8==0)
		{
		  
		  indix=0;
		  indix2++;
		  char* canvasname = new char[2];
		  sprintf(canvasname,"c%1d",indix2);
		  
		  fCanvasList->Add(new TCanvas(canvasname,"Traces ",600,800));
		  delete canvasname;
		  fCanvas=(TCanvas *) fCanvasList->Last();
		  fCanvas->Divide(2,4);
		  cout<<"new canvas"<<endl;
		  
		} 
	      
	      indix++;
	      
	      
	      
	      
	      char*hname = new char[15];
	      int idistance = (int)sta->fR_sf;
	      
	      int id = sta->fId;
	      
	      sprintf(hname,"PMT%1d_%04d_%04d",ipm+1,id,idistance);
	      
	      
	      cout<<hname<<endl;
	      histo=  new TH1F(hname,hname,nmax,0.,nmax);
	      // histo_mu=  new TH1F(hname,hname,nmax,0.,nmax);
	      delete hname;
	      
	      fHisto1->Add(histo);
	      // fHisto1->Add(histo_mu);
	      
	      
	      // adcsum=sta->fPMT[ipm].fADC;
	      for(int k=0; k<nmax; k++)
		if(sta->fPMT[ipm].fADC[0][k]!=0){
		  histo->Fill(25*k,sta->fPMT[ipm].fADC[0][k]/VEMPEAKVALUEINADC);
		  sum+=sta->fPMT[ipm].fADC[0][k];
	
		}
	      
	      //	      adcsum=sta->fPMT[ipm].fADC_mu;
	     //   for(int k=0; k<nmax; k++)
//  		if(sta->fPMT[ipm].fADC[0][k]!=0){
//  		  histo_mu->Fill(25*k,sta->fPMT[ipm].fADC[0][k]/VEMPEAKVALUEINADC);
//  		}
	      fCanvas->cd(indix);
	      histo->SetLineColor(kRed);
	      //   histo_mu->SetLineColor(kCyan);
	      histo->Draw();
	      //  histo_mu->Draw("SAME");
	      
	      cout<<"signal of station "<< sta->fId<<" = "<<sum/VEMCHARGEVALUEINADC<<" VEM"<<endl;
	      
	}
      
      
      if( indix%8==0)
	{
	  fCanvas->Print();
	}      
      
    }// end of loop on station
  
}






void Analyze::GetEvent(Int_t event)
{
  ofstream outfile;
  outfile.open("sum_cuve.out",ios::out);
  if (event>= fTree->GetEntries()) 
    {
      cout<<"Error in Analyze::GetShower() ; only " << (Int_t) fTree->GetEntries() << " events in Tree"<<endl;
      return;  
    }	
 
 
   cout<<"test"<<endl;
   fTree->GetEvent(event);
   
   fStList = fEvent->fHitStationList;
  
 
  
  
  cout << " Number of hit tanks ::" <<fStList.size()<< " "<<fEvent->fNombTank<<endl;
 
  vector<HitStation>::iterator sta;
  int i=0;
 int j=0;
 int k=0;
  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
    {
    if(sta->fT1Threshold ||sta->fT2ToT ||sta->fT2Threshold)
      // cout <<i<< "tank " <<sta->fId <<" pel = "<<sta->fNpe<<
// 	"\t nmu = "<<sta->fNmu<<" T1th "<<sta->fT1Threshold<<" ToT
// 	"<<sta->fT2ToT<<" T2th "<<sta->fT2Threshold<<endl;
    //  if(sta->fT2ToT ||sta->fT1Threshold)k++;
//      if(sta->fT2ToT ||sta->fT2Threshold)j++;
      outfile <<sta->fId <<"\t "<<sta->fR_sf<<"\t"<<sta->fNph<<"\t"<<sta->fNel<<"\t"<<sta->fNmu<<"\t"<<sta->fNpe_ph<<"\t"<<sta->fNpe_el<<"\t"<<sta->fNpe_mu<<" "<<(int)(sta->fNpe_ph/sta->fNph)<<" "<<(int)(sta->fNpe_el/sta->fNel)<<" "<<(int)(sta->fNpe_mu/sta->fNmu)<<" "<<sta->fNpe_mu/(sta->fNpe_el+sta->fNpe_ph)<<endl;

     cout <<i<< "\t tank " <<sta->fId <<"\t pel = "<<sta->fNpe<<         	"\t nmu = "<<sta->fNmu<<"\t nel = "<<sta->fNel<<"\t nph = "<<sta->fNph<<"\t dist = "<<sta->fR_sf<<"\t samp = "<<sta->fSampFact<<"\t time = "<<sta->fTime0<<"\t pemu= "<<sta->fNpe_mu<<"\t peel= "<<sta->fNpe_el<<"\t peph = "<<sta->fNpe_ph<<endl ; 
      i++;    
    }
 //  cout<<k<<"  stations passed a trigger "<<endl;
//   cout<<j<<"  stations passed a T2  "<<endl;
//  cout<<"Trigger T3  ";
//   if(fEvent->fT3Algo==1)cout<<"3TOT"<<endl;
//   else if(fEvent->fT3Algo==2)cout<<"4FOLD"<<endl;
//   else cout<< "FAILED"<<endl;

//   cout<<"Trigger T4  ";
//   if(fEvent->fT4Algo==1)cout<<"3TOT"<<endl;
//   else if(fEvent->fT4Algo==2)cout<<"4C1"<<endl;
//   else cout<< "FAILED"<<endl;
  
  
  int ista=0;
  double moymu=0;
  double moyel=0;
  double moyph=0;
  cout<<endl;
 //  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
//     {
//        if(sta->fNmu>0 && sta->fNel>8 && sta->fNph>100){
	
//  	cout <<sta->fId <<"  R= "<<sta->fR_sf<<
//  	  "  eff mu = "<<sta->fNpe_mu/(double)sta->fNmu<<
//  	  " nmu = "<<sta->fNmu<<" npemu = "<<sta->fNpe_mu<<
//  	  "  eff el = "<<sta->fNpe_el/(double)sta->fNel<<
//  	  " nel = "<<sta->fNel<<" npeel = "<<sta->fNpe_el<<
//  	  "  eff ph = "<<sta->fNpe_ph/(double)sta->fNph<<
//  	  " nph = "<<sta->fNph<<" npeph = "<<sta->fNpe_ph<<endl; 
//  	ista++; 
//  	moymu+=sta->fNpe_mu/(double)sta->fNmu;
//  	moyel+=sta->fNpe_el/(double)sta->fNel;
//  	moyph+=sta->fNpe_ph/(double)sta->fNph;
//        }
//     }
//   cout<<endl;
//   cout<<" average efficiency of particles "<<endl;
//   cout<<" mu= "<<moymu/(double)(ista)<<endl;
//   cout<<" el= "<<moyel/(double)(ista)<<endl;
//   cout<<" ph= "<<moyph/(double)(ista)<<endl;

 i=0;
  for(sta=fStList.begin();  sta!= fStList.end(); sta++)
    {
      //if(sta->fT2ToT ||sta->fT1Threshold)
      //if(((sta->fT2ToT==0) && (sta->fT1Threshold==0)) && ( sta->fR_sf < 2000 ))
      //{
      //cout <<i<< "\t tank " <<sta->fId <<"\t triggers = "<<sta->fT2ToT<<" "<<sta->fT2Threshold<<" \t time ="<<sta->fTime0<<"\t petot= "<<sta->fNpe<< "  dist_Sf= " << sta->fR_sf << endl ; 
      //if (sta->fR_sf < 300.)
      cout <<i<< "\t tank " <<sta->fId << "  dist_Sf= " << sta->fR_sf <<" \t time ="<<sta->fTime0<<"\t petot= "<<sta->fNpe<< endl ; 
      i++;
      //}    
    }



}
void Analyze::GetEvents()
{

  ofstream outfile;
  outfile.open("sum_cuve.out",ios::out);
  for (int event=0;event< fTree->GetEntries();event++) 
    {
      
      
      cout<<"test"<<endl;
      fTree->GetEvent(event);
      
      fStList = fEvent->fHitStationList;
      
      
      
      
      cout << " Number of hit tanks ::" <<fStList.size()<< " "<<fEvent->fNombTank<<endl;
      
      vector<HitStation>::iterator sta;
      int i=0;
      int j=0;
      int k=0;
      for(sta=fStList.begin();  sta!= fStList.end(); sta++){
      if(sta->fT1Threshold ||sta->fT2ToT ||sta->fT2Threshold)
	if(sta->fR_sf<2000)
	  outfile <<sta->fId <<"\t "<<sta->fR_sf<<"\t"<<sta->fNph<<"\t"<<sta->fNel<<"\t"<<sta->fNmu<<"\t"<<sta->fNpe_ph<<"\t"<<sta->fNpe_el<<"\t"<<sta->fNpe_mu<<" "<<(sta->fNpe_ph/sta->fNph)<<" "<<(sta->fNpe_el/sta->fNel)<<" "<<(sta->fNpe_mu/sta->fNmu)<<" "<<sta->fNpe_mu/(sta->fNpe_el+sta->fNpe_ph)<<endl;
      
    }
      

    }

}

void Analyze::GetEvents(int ista)
{

  ofstream outfile;
  int nmu,nel,nph,npemu,npeph,npeel;
  float rael,ramu,raph,ramuel,nev;
  outfile.open("sum_cuve.out",ios::out);
  nmu=0;
  nel=0;
  nph=0;
  npeph=0;
  npeel=0;
  npemu=0;
  rael=0;
  ramu=0;
  raph=0;
  ramuel=0;
  nev=(float)fTree->GetEntries();
  for (int event=0;event< fTree->GetEntries();event++) 
    {
      
      
      cout<<"test"<<endl;
      fTree->GetEvent(event);
      
      fStList = fEvent->fHitStationList;
      
      
      
      
      cout << " Number of hit tanks ::" <<fStList.size()<< " "<<fEvent->fNombTank<<endl;
      
      vector<HitStation>::iterator sta;
      int i=0;
      int j=0;
      int k=0;
      for(sta=fStList.begin();  sta!= fStList.end(); sta++){
	if(sta->fId!=ista)continue;
      if(sta->fT1Threshold ||sta->fT2ToT ||sta->fT2Threshold)
	if(sta->fR_sf<2000)
	  //  outfile <<sta->fId <<"\t "<<sta->fR_sf<<"\t"<<sta->fNph<<"\t"<<sta->fNel<<"\t"<<sta->fNmu<<"\t"<<sta->fNpe_ph<<"\t"<<sta->fNpe_el<<"\t"<<sta->fNpe_mu<<" "<<(sta->fNpe_ph/sta->fNph)<<" "<<(sta->fNpe_el/sta->fNel)<<" "<<(sta->fNpe_mu/sta->fNmu)<<" "<<sta->fNpe_mu/(sta->fNpe_el+sta->fNpe_ph)<<endl;
      nmu+=sta->fNmu;
      nel+=sta->fNel;
      nph+=sta->fNph;
      npeph+=sta->fNpe_ph;
      npemu+=sta->fNpe_mu;
      npeel+=sta->fNpe_el;
      rael+=sta->fNpe_el/sta->fNel;
      ramu+=sta->fNpe_mu/sta->fNmu;
      raph+=sta->fNpe_ph/sta->fNph;
      ramuel+=sta->fNpe_mu/(sta->fNpe_el+sta->fNpe_ph);
      

    }
      

    }
  outfile<<ista<<" "<<nph/nev<<" "<<nel/nev<<" "<<nmu/nev<<" "<<npeph/nev<<" "<<npeel/nev<<" "<<npemu/nev<<" "<<raph/nev<<" "<<rael/nev<<" "<<ramu/nev<<" "<<ramuel/nev<<endl;
}



void Analyze::PlotDistances()
{
  TCanvas *CanLTP= new TCanvas("Distances","Distances",600,800);
  CanLTP->SetFillColor(0);
  
  vector<HitStation>::iterator sta;  
  TH1F* dist=new TH1F("Distances","Distances in Shower Frame",80,0.,4000); 
  
  
  for (int i=0;i< fChain->GetEntries() ; i++) 
    {
      fChain->GetEntry(i);
      fStList = fEvent->fHitStationList;
      
      for(sta=fStList.begin();  sta!= fStList.end(); sta++)
	{
	  dist->Fill(sta->fR_sf,1.);
	  if (sta->fR_sf<200.)
	    cout <<"tank " <<sta->fId << "  dist_Sf= " << sta->fR_sf << endl ; 
	}
    }
  dist->Draw();

}




void Analyze::PlotPMTBalance()
{
  vector<HitStation>::iterator sta; 
  if(fChain!=0)
    { 
      delete fChain;
      fChain=0;
    }
  fChain= new TChain("TreeEvent");
  
  
  ifstream in1;
  string Sim="/windows/partage/Showers";
  //  string  Primary="proton/";
  string  Primary="Olivier/energy1";
  //string  Theta="0.00d/";
  const int ntheta=1;
  const int nenergy=1;
  const int nphi=1;
 //   string Theta[ntheta]={ "0.00","25.84","36.87","45.57","53.13","60.00"};
//    string Energy[nenergy]={"2.512e+08","6.31e+08","1e+09","1.585e+09","2.512e+09",
//  		     "3.981e+09","6.31e+09","1.585e+10","3.981e+10","6.31e+10",
//  		     "1e+11","1.585e+11"};
//    int Phi[nphi]={54,126,198,270,342};



  int Phi=126;
  string Theta="56.63";
  string Energy="1.585e+09";
  string DebutFichier="Sim_DAT00*";
  string tampon="tampon.dat";  
  
  string commandrm="rm tampon.dat";
  char * readfilename1= new char[250];
  char * readfilename2= new char[250];; 
  char * commandls= new char[250];;
  
  for(int i=0;i<ntheta;i++)
    for(int j=0;j<nenergy;j++)
      for(int k=0;k<nphi;k++){
	
	
	//sprintf(commandls,"ls %s/%s/%sd/%s/%dd/%s > %s",Sim.c_str(),Primary.c_str(),Theta[i].c_str(),Energy[j].c_str(),Phi[k],DebutFichier.c_str(),tampon.c_str());
	sprintf(commandls,"ls %s/%s/%s/%s/%dd/%s > %s",Sim.c_str(),Primary.c_str(),Theta.c_str(),Energy.c_str(),Phi,DebutFichier.c_str(),tampon.c_str());
	cout << commandls <<endl;
	system(commandls);
	
	in1.open(tampon.c_str(),ios::in);
	
	in1 >> readfilename1;
	cout << readfilename1 << endl;
	fChain->Add(readfilename1);
        in1 >> readfilename2;
	cout << readfilename2 << endl;
	fChain->Add(readfilename2);
	
	in1.close();
	
	cout << commandrm << endl;
	system(commandrm.c_str());
	
      }
  
  
  
  delete [] commandls;
  cout << "toto" << endl;
  delete [] readfilename1;
  cout << "toto2" << endl;
  delete [] readfilename2;
  cout << "toto3" << endl;
  //delete in1;
  cout << "toto4" << endl;
  if(fBranch!=0)
    { 
      delete fBranch;
      fBranch=0;
    }
  
  fBranch= fChain->GetBranch("Event");
  fChain->SetBranchAddress("event",&fEvent);
  
  cout << " \t * There is " << fChain->GetEntries() << " events in the input file "<<endl;
  
  fChain->GetEntry(0);
  fChain->Print();
 
 TFile * file=new TFile("testsignal2.root","RECREATE");
 TNtuple * nt=new TNtuple("testsignal","testsignal","id:E:theta:phi:ntanks:tankid:pm1:pm2:pm3");


  for(int i=0;i<fChain->GetEntries();i++){
    
    fChain->GetEntry(i);
    fStList = fEvent->fHitStationList;
    
    for(sta=fStList.begin();  sta!= fStList.end(); sta++)
      {
	float pm1=0;
	float pm2=0;
	float pm3=0;

	for(int j=0;j<MAXNUMBEROFADCBINS;j++){
	  pm1+=sta->fPMT[0].fADC[0][j];
	  pm2+=sta->fPMT[1].fADC[0][j];
	  pm3+=sta->fPMT[2].fADC[0][j];
	}
	nt->Fill((float)i,(float)fEvent->fEnergy,(float)fEvent->fTheta,(float)fEvent->fAzim,(float)fEvent->fNombTank,(float)sta->fId,pm1/177.2,pm2/177.2,pm3/177.2);
	
      }
  }
  
   file->Write();
   file->Close();
   
   return;
  
  
}


void Analyze::PrintTree()
{  
  fTree->Print();
}



void Analyze::ReadFile(char * filename)
{
  
  if(fTree!=0)
    { 
      delete fTree;
      fTree=0;
    }
  
  fFile = new TFile(filename,"READ");	// Open file	
  fFile->ls();
  
  fTree = (TTree *) fFile->Get("TreeEvent");      // Get Tree
  
  fTree->Print(); 
  
  if(fBranch!=0)
    { 
      delete fBranch;
      fBranch=0;
    }
  
  
  
  
  
  fBranch= fTree->GetBranch("event");
  fTree->SetBranchAddress("event",&fEvent);
  
  cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
  
  // fTree->GetEntry(0);
  
  return;
}





//--------------------------------------------------------------------









void Analyze::TestAnodeDynode()
{

  map<int,double> pm_hi;
  map<int,double> pm_lo;
  ofstream outfilehi,outfilelo;
  HitStation* sta;
  int nbsta;
 
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  int nev=(int)fTree->GetEntries();
  cout<<nev<<" events to use"<<endl;
 
  outfilehi.open("pmhi.dat",ios::out);
  outfilelo.open("pmlo.dat",ios::out);

  for(int ne=0; ne<=nev; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      
      nbsta= fStList.size();
      for(int ista=0;ista<nbsta;ista++){
	sta=&(fStList[ista]);
	
	
	if(sta->fNpe>2000) 
	  {
	    pm_hi =sta->fPMT[0].fPMTSignal_hi;
	    pm_lo =sta->fPMT[0].fPMTSignal_lo;
	    double samp=sta->fSampFact;
	    int firstbin=-100;
	   
	    for(map<int,double>::const_iterator slot=pm_hi.begin();
		slot!=pm_hi.end(); slot++)
	      {
		if(((*slot).first<10000))
		  {
		    firstbin=(*slot).first;
		    break;
		  }
	      }
		
	    if(firstbin>0)
	      for(int k=firstbin; k<5000+firstbin; k++){
		if(pm_hi.count(k))
		 outfilehi<<" "<<pm_hi[k]/samp<<" "<<pm_hi[k]/samp<<" "<<pm_hi[k]/samp<<endl;
		else
		  outfilehi<<" 0 0 0"<<endl;
	      }
	    
	     firstbin=-100;
	     for(map<int,double>::const_iterator slot=pm_lo.begin();
		 slot!=pm_lo.end(); slot++)
	       {
		 if(((*slot).first<10000))
		   {
		     firstbin=(*slot).first;
		     break;
		   }
	       }
	     if(firstbin>0)
	       for(int k=firstbin; k<5000+firstbin; k++){
		 if(pm_lo.count(k))
		   outfilelo<<" "<<pm_lo[k]/samp<<" "<<pm_lo[k]/samp<<" "<<pm_lo[k]/samp<<endl;
		 else
		   outfilelo<<" 0 0 0"<<endl;
	       }
	  }
	
	
      }
    }
  outfilelo.close();
  outfilehi.close();
}
	




void Analyze::TestNewTOT()
{
  HitStation* sta;
  char *filename= new char[40];
  float LowThreshold = 0.2*24.;
  float UpThreshold = 0.4*24;
  int WindowSize = 120;
  int NumberOfRealisationL = 10;
  int nentry = (int)fTree->GetEntries();
  int n;
  ofstream out1;
  
  out1.open("totut.dat",ios::out);
  
  
  
  
  for(int upperadd=0; upperadd<2; upperadd++){
    UpThreshold+=0.2*24;
    if( upperadd ==4 )UpThreshold=1000;
    for(int realadd=0;realadd<2;realadd++)
      {
	int NumberOfRealisation=NumberOfRealisationL+realadd*2;
	n=0;
	for(int i = 0; i < nentry; i++){
	  
	  
	  fTree->GetEntry(i);
	  fStList = fEvent->fHitStationList;
	  sta=&(fStList[0]);
	  if(i%1000==0) cout<< i << " particles "<<endl;
	  
	 
	  
	  
	  int count=0;
	  
	  for(int ipm=0; ipm<NPM; ipm++)
	    {
	      
	      
	      for(int ibin=0;ibin<400;ibin++){ 
		int count2 = 0;
		int begin = 0;
		if(sta->fPMT[ipm].fADC[0][ibin]>LowThreshold && sta->fPMT[ipm].fADC[0][ibin]<UpThreshold){
		  begin = ibin;
		  
		  
		  do{
		    if(sta->fPMT[ipm].fADC[0][ibin] > LowThreshold && sta->fPMT[ipm].fADC[0][ibin]<UpThreshold) count2+=1;
		    ibin++;
		  }while(ibin<(begin+WindowSize) && ibin<400);
		}
		
		
		
		if(count2>=NumberOfRealisation)
		  {
		    if(upperadd!=4) out1<<" TOT sur PM#"<<ipm+1<< "  muon "<<i<<endl;
		    // out1<<endl;
		    count+=1;
		    break;
		  }
	      }
	      
	      
	    }
	  if(count>=2)
	    {
	      n++;
	      out1<<"TOT SUR DOUBLE MU Numero = "<<i<<" "<<n<<endl;
	    }
	 
	}
	out1<<WindowSize<<" "<<UpThreshold<<" "<<NumberOfRealisation<< " "<<n<<endl;
      }
  }
  
  
  out1.close();
  delete filename;
}

/*
void Analyze::DrawMuons()
{
  map<int,double>* timeprof,*timeprofsum;
  map<int,double>* pmtsignalsum;
  double fefiltersum[MAXNUMBEROFTIMEBINS];
  short adcsum[2][MAXNUMBEROFADCBINS];
  int nmu,nper;
  double npesum,npe[NPM],npemax;
  HitStation* sta;
  TProfile* profpe[NPM],*profpesum,*profpmsum,*proffesum,*profadcsum;
  TH1F* hpe[NPM],*hpesumch,*hpmsumch,*hfesumch,*hadcsumch;
  TH1F*hpesumpk,*hpmsumpk,*hfesumpk,*hadcsumpk;
  TH1F*hpesumatopk,*hpmsumatopk,*hfesumatopk,*hadcsumatopk;
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);

 
  hpesumch =  new TH1F("Pelsumcharge","Pelsumcharge",600,0.,600);
  hpesumpk =  new TH1F("Pelsumpeak","Pelsumpeak",100,0.,100);
  hpesumatopk =  new TH1F("Pelsumatopk","Pelsumatopk",100,0.,100);
  profpesum = new TProfile("profpelsum","profpelsum",500,0,500,0,500);
  hpmsumch =  new TH1F("Pmsumcharge","pmsumcharge",500,0.,50);
  hpmsumpk =  new TH1F("Pmsumpeak","pmsumpeak",100,0.,1);
  hpmsumatopk =  new TH1F("Pmsumatopk","pmsumatopk",400,0.,200);
  profpmsum = new TProfile("profpmsum","profpmsum",500,0,500,0,500);
  hfesumch =  new TH1F("fesumcharge","fesumcharge",500,0.,50);
  hfesumpk =  new TH1F("fesumpeak","fesumpeak",100,0.,1);
  hfesumatopk =  new TH1F("fesumatopk","fesumatopk",400,0.,200);
  proffesum = new TProfile("proffesum","proffesum",500,0,500,0,500);
  hadcsumch =  new TH1F("adcsumcharge","adcsumcharge",600,0.,600);
  hadcsumpk =  new TH1F("adcsumpeak","adcsumpeak",600,0.,600);
  hadcsumatopk =  new TH1F("adcsumatopk","adcsumatopk",100,0.,10);
  profadcsum = new TProfile("profadcsum","profadcsum",500,0,50,0,50);

  for(int ipm=0;ipm<NPM;ipm++){
    char*hname = new char[15];
    sprintf(hname,"Photoelpm%d",ipm+1);
    hpe[ipm]=  new TH1F(hname,hname,200,0.,200);
    profpe[ipm] = new TProfile(hname,hname,500,0,500,0,500);
    delete hname;
  }

  nmu=(int)fTree->GetEntries();
  cout<<nmu<<" muons to analyze"<<endl;



  for(int ne=0; ne<nmu; ne++)
    {
      fTree->GetEvent(ne);
      if(ne%1000==0) cout<< ne << " particles "<<endl;
      fStList = fEvent->fHitStationList;
      if(fStList.size()!=1) cout<<" Error nb of tanks "<<fStList.size()<<endl;
      sta=&(fStList[0]);
  

      // work on sum traces  Photoelectrons
      npesum=0;
      npemax=0;
      timeprofsum =&(sta->fTimeProfile);
      nper=sta->fNpe;
      // loop on timebins
     
      for(int k=0; k<1000; k++){
	if(timeprofsum->count(k)){
	  profpesum->Fill(k,(*timeprofsum)[k]);
	  npesum+=(*timeprofsum)[k] ;
	  if((*timeprofsum)[k]>npemax) npemax=(*timeprofsum)[k];
	} 
	else
	  profpesum->Fill(k,0);
      }// end of loop on k bins
      hpesumch->Fill(npesum);
      hpesumpk->Fill(npemax);
      if(npemax!=0)hpesumatopk->Fill(npesum/npemax);
      
      // work with pmtsignal  
      npesum=0;
      npemax=0;
      pmtsignalsum =&(sta->fPMTSignal_hi);
      for(int k=0; k<1000; k++){
	  if(pmtsignalsum->count(k)){
	    profpmsum->Fill(k,(*pmtsignalsum)[k]);
	    npesum+=(*pmtsignalsum)[k] ;
	    if((*pmtsignalsum)[k]>npemax) npemax=(*pmtsignalsum)[k];
	  } 
	  else
	    profpmsum->Fill(k,0);
      }// end of loop on k bins
      hpmsumch->Fill(npesum);
      hpmsumpk->Fill(npemax);
      if(npemax!=0)hpmsumatopk->Fill(npesum/npemax);

      // work with fesignal  
      npesum=0;
      npemax=0;
      fefiltersum =sta->fFEFilter_hi;
      
      // new loop on timebins
      for(int k=0; k<1000; k++){
	if(fefiltersum[k]!=0){
	  proffesum->Fill(k,fefiltersum[k]);
	  if(fefiltersum[k]>npemax) npemax=fefiltersum[k];
	  npesum+=fefiltersum[k] ;
	} 
	else
	  proffesum->Fill(k,0);
      }// end of loop on k bins
      hfesumch->Fill(npesum);
      hfesumpk->Fill(npemax);
      if(npemax!=0)hfesumatopk->Fill(npesum/npemax);

      // work with adcsignal  
      npesum=0;
      npemax=0;
      adcsum =sta->fADC;
      
      // new loop on timebins
      for(int k=0; k<1000; k++){
	if(adcsum[0][k]!=0){
	  profadcsum->Fill(k,adcsum[0][k]);
	  npesum+=adcsum[0][k] ;
	  if(adcsum[0][k]>npemax) npemax=adcsum[0][k];
	} 
	else
	  profadcsum->Fill(k,0);
     }// end of loop on k bins
      hadcsumch->Fill(npesum);
      hadcsumpk->Fill(npemax);
      if(npemax!=0)hadcsumatopk->Fill(npesum/npemax);
      
      // work on individual PMT traces
      for(int ipm=0; ipm<NPM;ipm++){
	npe[ipm]=0;
	timeprof =&(sta->fPMT[ipm].fTimeProfile);
	nper=sta->fPMT[ipm].fNpe;
	// new loop on timebins
	for(int k=0; k<1000; k++){
	  if(timeprof->count(k)){
	    profpe[ipm]->Fill(k,(*timeprof)[k]);
	    npe[ipm]+=(*timeprof)[k] ;
	  } 
	  else
	    profpe[ipm]->Fill(k,0);
	}// end of loop on k bins
	hpe[ipm]->Fill(npe[ipm]);
	
      }//end of loop on PM
     
    }
  fCanvasList->Add(new TCanvas("photoel"," photoel ",800,800));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,3);
 
  for (int ipm=0;ipm<NPM;ipm++){
    fCanvas->cd(2*ipm+1);
    hpe[ipm]->SetLineColor(kRed);
    hpe[ipm]->Draw(); 
    fCanvas->cd(2*ipm+2);
    profpe[ipm]->SetLineColor(kRed);
    profpe[ipm]->Draw();  
  }

  fCanvasList->Add(new TCanvas("profiles"," profiles ",800,800));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,4);
  fCanvas->cd(1);
  hpesumch->SetLineColor(kRed);
  hpesumch->Draw();
  fCanvas->cd(2);
  profpesum->SetLineColor(kRed);
  profpesum->Draw();
  fCanvas->cd(3);
  hpmsumch->SetLineColor(kRed);
  hpmsumch->Draw();
  fCanvas->cd(4);
  profpmsum->SetLineColor(kRed);
  profpmsum->Draw();
  fCanvas->cd(5);
  hfesumch->SetLineColor(kRed);
  hfesumch->Draw();
  fCanvas->cd(6);
  proffesum->SetLineColor(kRed);
  proffesum->Draw();
  fCanvas->cd(7);
  hadcsumch->SetLineColor(kRed);
  hadcsumch->Draw();
  fCanvas->cd(8);
  profadcsum->SetLineColor(kRed);
  profadcsum->Draw();

 fCanvasList->Add(new TCanvas("muAnalyse"," mu analyse ",800,800));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,4);
  fCanvas->cd(1);
  hpesumpk->SetLineColor(kRed);
  hpesumpk->Draw();
  fCanvas->cd(2);
  hpesumatopk->SetLineColor(kRed);
  hpesumatopk->Draw();
  fCanvas->cd(3);
  hpmsumpk->SetLineColor(kRed);
  hpmsumpk->Draw();
  fCanvas->cd(4);
  hpmsumatopk->SetLineColor(kRed);
  hpmsumatopk->Draw();
  fCanvas->cd(5);
  hfesumpk->SetLineColor(kRed);
  hfesumpk->Draw();
  fCanvas->cd(6);
  hfesumatopk->SetLineColor(kRed);
  hfesumatopk->Draw();
  fCanvas->cd(7);
  hadcsumpk->SetLineColor(kRed);
  hadcsumpk->Draw();
  fCanvas->cd(8);
  hadcsumatopk->SetLineColor(kRed);
  hadcsumatopk->Draw();
}}
*/
/*

void Analyze::Help()
{
  cout << endl;
  cout << "Available commands: " << endl;
  cout << endl;
  cout <<"Readfile('foo.root')	    : load  output file  "<<endl;
  cout <<"GetShower(n)		    : Load nth event in memory" << endl;
  cout <<"GetAllShowers		    : Load all event in memory" << endl;
  cout <<"PrintTree()		    : fTree->Print()" << endl;
  cout <<"DrawTraces(ipm,icomp)"<<endl;
  cout <<"DrawOneTrace(nutank in vect,ipm,icomp)" <<endl;
  cout <<"DrawSumTraces()"<<endl;
  cout <<"StatShowers()"<<endl;
  cout<<"Shape()                    : plot and compare identification parameters"<<endl;
 
}
void Analyze::DrawEvent()
{
  ifstream in1;
  in1.open("station.dat",ios::in);
  int numb = fTanks.size();
  int n = 0;
  int numb2 = 40 - numb;


  vector<Tank>::iterator sta;
  fCanvasList->Add(new TCanvas("Event","Event Map",600,800));
  fCanvas=(TCanvas *) fCanvasList->Last();



  fCanvas->Clear();
  fCanvas->Range(-6000,-6000,6000,8000);
  fCanvas->SetBorderSize(2);
  fCanvas->cd();
  TEllipse *nonstation[numb2];
  //TText *nonstationid[numb2];
  TEllipse *station[numb];
  TText *stationid[numb];
  TMarker *Core;
  Core = new TMarker();
  Core->SetMarkerSize(1);
  Core->SetMarkerColor(2);
  Core->SetMarkerStyle(29);
  Core->SetX(fShower->fxcenter);
  Core->SetY(fShower->fycenter);
  cout<<"CACA "<<fShower->fxcenter<<" " <<fShower->fycenter<<" "<<fShower->fangle<<" "<<fShower->fphi<<endl;
  Core->Draw();
  fCanvas->Modified();
  fCanvas->Update();
  cout<<"merde "<<endl;
  TLine *ShowerAxis = new TLine();
  ShowerAxis->SetLineColor(2);
  ShowerAxis->SetLineWidth(2);
  ShowerAxis->SetX1(fShower->fxcenter);
  ShowerAxis->SetX2(fShower->fxcenter+2500*cos(fShower->fphi)*sin(fShower->fangle/180*3.1416));
  ShowerAxis->SetY1(fShower->fycenter);
  ShowerAxis->SetY2(fShower->fycenter+2500*sin(fShower->fphi)*sin(fShower->fangle/180*3.1416));
  cout<<fShower->fangle<<" "<<fShower->fphi * 180/3.1416<<endl;
  for(sta=fTanks.begin();  sta!= fTanks.end(); sta++)
    {
      char* name = new char[5];
      sprintf(name,"%1d",sta->fid);
      double size = sqrt(sta->fphotoel_tot[NPM]/4)*5;
      if (size<200.)
        size=200.;
      if(sta->fid!=49 && sta->fid!=64)
        {
	 
          station[n]=new TEllipse(sta->fxstat,sta->fystat,size,size,0,360,0);
          stationid[n]=new TText(sta->fxstat,sta->fystat,name);
        }
      else
        {
          if(sta->fid==49)
            {
              station[n]=new TEllipse(sta->fxstat,sta->fystat,size,size,0,180,90);
              stationid[n]=new TText(sta->fxstat-250,sta->fystat,name);
            }
          else
            {
              station[n]=new TEllipse(sta->fxstat,sta->fystat,size,size,0,180,270);
              stationid[n]=new TText(sta->fxstat+250,sta->fystat,name);
            }

        }
      delete name;
      stationid[n]->SetTextAlign(22);
      //     stationid[n]->SetTextColor(9);
      stationid[n]->SetTextSize(0.024);
      stationid[n]->SetTextAngle(40);
      if(sta->fflagupdown==1)station[n]->SetLineColor(kRed);
      if(sta->fflagupdown==-1)station[n]->SetLineColor(kBlue);
      if(sta->fflagupdown==0)station[n]->SetLineColor(kMagenta);
      if(sta->fflagupdown==1)stationid[n]->SetTextColor(kRed);
      if(sta->fflagupdown==-1)stationid[n]->SetTextColor(kBlue);
      if(sta->fflagupdown==0)stationid[n]->SetTextColor(kMagenta);

      //station[n]->SetLineColor(2);
      station[n]->SetLineWidth(2);
      station[n]->SetLineStyle(2);
      station[n]->Draw();
      stationid[n]->Draw();
      n++;


    }

  int id;
  double x, y;
  int p = 0;
  int is;
  ShowerAxis->Draw();
  do
    {

      in1 >> id >> y>> x;
      if (in1.eof()) break;
      is = 0;
      for(sta=fTanks.begin();  sta!= fTanks.end(); sta++)
        {
          if(sta->fid==id)
            is=1;
        }
      if(is==0)
        {
          //  char* name2 = new char[5];
          // sprintf(name2,"%1d",id);

          if(id!=49 && id!=64)
            {
              nonstation[p] = new TEllipse(x-459630.,y-6.08276E+06,200,200,0,360,0);
	     

              //  nonstationid[p] = new TText(x-459630.,y-6.08276E+06,name2);

            }
          else
            {
              if(id==49)
                {
                  nonstation[p] = new TEllipse(x-459630.,y-6.08276E+06,200,200,0,180,90);

                  //  nonstationid[p] = new TText(x-459630.-500,y-6.08276E+06,name2);
                }
              else
                {
                  nonstation[p] = new TEllipse(x-459630.,y-6.08276E+06,200,200,0,180,270);
                  // nonstationid[p] = new TText(x-459630.+500,y-6.08276E+06,name2);
                }
            }

          //  delete name2;
          // nonstationid[p]->SetTextAlign(22);
          //nonstationid[p]->SetTextSize(0.03);
          // nonstationid[p]->SetTextAngle(40);
          nonstation[p]->SetLineWidth(2);
          nonstation[p]->SetLineStyle(2);
          nonstation[p]->Draw();
          //  nonstationid[p]->Draw();

        }
    }while (!in1.eof());
  cout<<" canvas dessine"<<endl;
 
  TLatex *L = new TLatex();
  cout<<"coucou1"<<endl;
  L->SetTextAlign(12);
  cout<<"coucou1"<<endl;
  L->SetTextSize(0.04);
  cout<<"coucou1"<<endl;
  char* value = new char[40];
  sprintf(value,"E = %g 10^{18} eV",fShower->fenergy);
  L->DrawLatex(-5000,7500,value);
  delete value;
  cout<<"coucou1"<<endl;
  char* value2 = new char[40];
  sprintf(value2,"#theta = %g#circ",fShower->fangle);
  L->DrawLatex(-5000,6900,value2);
  delete value2;
  cout<<"coucou1"<<endl;
  char* value3 = new char[40];
  sprintf(value3,"#phi = %g#circ",fShower->fphi*180/3.14159);
  L->DrawLatex(-5000,6300,value3);
  delete value3;
  in1.close();
 
}
*/

