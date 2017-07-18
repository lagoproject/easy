#include "AnalyzeStar.h"
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
using namespace std;

// NOTE: this code is not used in LAGO

ClassImp(AnalyzeStar)


  AnalyzeStar::AnalyzeStar()
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

AnalyzeStar::~AnalyzeStar()
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

void AnalyzeStar::CompareProfiles()
{
  ifstream files;
  char * filename;
  TFile * file;
  filename = new char[40];
  int nsta;
// short adcsum[2][MAXNUMBEROFADCBINS];
// short adcsum_ns[MAXNUMBEROFADCBINS];
  int nmax = MAXNUMBEROFADCBINS;
  int sumtot,sumsat;
  TProfile* lat[2];
  TProfile* latsat[2];
  fCanvasList->Add(new TCanvas("comp","comp sampling ",800,800));
  fCanvas=(TCanvas *) fCanvasList->Last();
  fCanvas->Divide(2,2);

 
  latsat[0] = new TProfile("latsat1","latsat1",2500,0,2500,0,100000); 
  lat[0] = new TProfile("lat1","lat1",2500,0,2500,0,100000); 
  latsat[1] = new TProfile("latsat2","latsat2",2500,0,2500,0,100000); 
   lat[1] = new TProfile("lat2","lat2",2500,0,2500,0,500000); 
  TProfile* dif = new TProfile("dif","dif",1500,0,1500,-10000,10000); 
  TProfile* difsat = new TProfile("difsat","difsat",1500,0,1500,-10000,10000); 
   TProfile* samp = new TProfile("samp","samp",1500,0,1500,0,10000); 
    TProfile* npart = new TProfile("npart","npart",1500,0,1500,0,100000); 
   TProfile* nmu = new TProfile("nmu","nmu",1500,0,1500,0,100000); 
   TProfile* nel = new TProfile("nel","nel",1500,0,1500,0,100000); 
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetLineColor(kRed);
  cout<<"lecrures des ficheirs "<<endl;
 
  files.open("compare.inp",ios::in);
 
  for (int i=0;i<2;i++)
    {
      files>>filename;
      cout<<filename<<endl;
      file = 0;
      fTree=0;
      fBranch=0;
      fEvent=0;
      
      
      file = new TFile(filename,"READ");	 // Open file	
      file->ls();
     
      fTree = (TTree *) file->Get("TreeEvent");      // Get Tree
      
      fTree->Print(); 
      
       
      fBranch= fTree->GetBranch("HitStation");
      fTree->SetBranchAddress("hstation",&fHitStation);
  

     
      cout << " \t * There is " << fTree->GetEntries() << " stations in the input file "<<endl;
      
    
      nsta=(int)fTree->GetEntries();
      
 
      
      for(int ns=0; ns<nsta; ns++)
	{
	  if(i==0)
	    {
	      samp->Fill(fHitStation->fR_sf,fHitStation->fSampFact);
	      npart->Fill(fHitStation->fR_sf,fHitStation->fNmu
			  +fHitStation->fNel +fHitStation->fNph);
	      nmu->Fill(fHitStation->fR_sf,fHitStation->fNmu);
	      nel->Fill(fHitStation->fR_sf,fHitStation->fNmu+
			fHitStation->fNel +fHitStation->fNph);

	    }
	  
	  sumtot=0;
	  sumsat=0;
// 	  	  for(int ii = 0;ii<MAXNUMBEROFADCBINS;ii++)
// 	    {
// 	      adcsum_ns[ii]=0;
// 	      for(int j = 0;j<2;j++)
// 		{
// 		  adcsum[j][ii]=0;
		  
// 		}
// 		}
	  
	  
	  cout<<"avant event"<<endl;
	  fTree->GetEvent(ns);
	
	  cout << "\t tank " <<fHitStation->fId <<"\t pel = "<<fHitStation->fNpe<<
	    "\t nmu = "<<fHitStation->fNmu<<"\t dist = "<<fHitStation->fR_sf<<endl ; 
	  
	
	  
//adcsum=fHitStation->fADC;
	  for(int k=0; k<nmax; k++)
	    sumsat+=(int)fHitStation->fADC[k];
	  latsat[i]->Fill(fHitStation->fR_sf,sumsat);
	  
// adcsum_ns=fHitStation->fADC_ns;
	  for(int k=0; k<nmax; k++)
	    sumtot+=fHitStation->fADC_ns[k];
	  lat[i]->Fill(fHitStation->fR_sf,sumtot);
	  cout<<sumsat<<" "<<sumtot<<endl;
	  
	}
	  
      fCanvas->cd(1);
      lat[i]->SetMarkerColor(kRed);
      if(i==0)lat[i]->Draw();
      else lat[i]->Draw("same");
	  
      fCanvas->cd(2);
      latsat[i]->SetMarkerColor(kCyan);
      
      
      if(i==0)latsat[i]->Draw();
      else latsat[i]->Draw("SAME");
      
      file->Close(0) ;
    }//end of loop on files
  
  for(int k=300; k<2000; k++)
   {
     double  lat0=lat[0]->GetBinContent(k);
     double  lat1=lat[1]->GetBinContent(k);
     double  latsat0=latsat[0]->GetBinContent(k);
     double  latsat1=latsat[1]->GetBinContent(k);
     
     if(lat0>0)
       {
       dif->Fill(k,(lat0-lat1)/lat0);
       cout<<k<<" "<<lat0<<" "<<lat1<<" "<<(lat0-lat1)/lat0<<endl;
       
       }
     if(latsat0>0)
       {
	 difsat->Fill(k,(latsat0-latsat1)/lat0);
	 cout<<k<<" "<<latsat0<<" "<<latsat1<<" "<<(latsat0-latsat1)/latsat0<<endl;
       }
   }

 fCanvas->cd(3);
 dif->Draw();
 fCanvas->cd(4);
 difsat->Draw();

 fCanvasList->Add(new TCanvas("samp"," sampling ",800,800));
 fCanvas=(TCanvas *) fCanvasList->Last();
 fCanvas->Divide(2,2);
 fCanvas->cd(1);
 samp->Draw();
 fCanvas->cd(2);
 npart->Draw();
 fCanvas->cd(3);
 nmu->Draw();
 fCanvas->cd(4);
 nel->Draw();

}


void AnalyzeStar::ComputeNbotDistr(int energy,int tta)
{
  double dist1=100.;
  double dist2=dist1+100.;
  double dist3=dist2+100.;
  double dist4=dist3+100.;
  double delta=50.;
  double mean1,mean2,mean3,mean4;
  double rms1,rms2,rms3,rms4;

  TH1D *Dist1; 
  TH1D *Dist2; 
  TH1D *Dist3; 
  TH1D *Dist4;
  
  char * outputfilename = new char[120];  
  sprintf(outputfilename,"condladibot%d_%d.txt",energy,tta); 
  ofstream out1;
  out1.open(outputfilename,ios::out);
  
  for (int j=0;j<5;j++)
    {
      Dist1 = new TH1D("Dist1","d=100m",100,0.,60);
      Dist2 = new TH1D("Dist2","d=200m",100,0.,60);
      Dist3 = new TH1D("Dist3","d=300m",100,0.,60);
      Dist4 = new TH1D("Dist4","d=400m",100,0.,60);
      dist1=100.+400.*j;
      dist2=dist1+100.;
      dist3=dist2+100.;
      dist4=dist3+100.;
  
      for (int i=0;i< fChain->GetEntries() ; i++) 
	{
	  fChain->GetEntry(i);
	  
	  
	  if (fHitStation->fT2ToT)
	    {
	      if ( (fHitStation->fR_sf > (dist1-delta)) && (fHitStation->fR_sf< (dist1+delta)) )
		{
		  Dist1->Fill((double)(fHitStation->fNbOfBinsOverThres));
		}
	      
	      if ((fHitStation->fR_sf> (dist2-delta)) && (fHitStation->fR_sf< (dist2+delta)))
		{
		  Dist2->Fill((double)(fHitStation->fNbOfBinsOverThres));
		}
	      
	      if ((fHitStation->fR_sf> (dist3-delta)) && (fHitStation->fR_sf< (dist3+delta)))
		{
		  Dist3->Fill((double)(fHitStation->fNbOfBinsOverThres));
		} 
	      
	      if ((fHitStation->fR_sf> (dist4-delta)) && (fHitStation->fR_sf< (dist4+delta)))
		{
		  Dist4->Fill((double)(fHitStation->fNbOfBinsOverThres));
		} 
	    }
	}
      mean1=Dist1->GetMean();
      rms1=Dist1->GetRMS();
      mean2=Dist2->GetMean();
      rms2=Dist2->GetRMS();
      mean3=Dist3->GetMean();
      rms3=Dist3->GetRMS();
      mean4=Dist4->GetMean();
      rms4=Dist4->GetRMS();	
      
      
    }
  
  out1 << mean1 << "  " << rms1 << endl;
  out1 << mean2 << "  " << rms2 << endl;	
  out1 << mean3 << "  " << rms3 << endl;	
  out1 << mean4 << "  " << rms4 << endl;		

  delete Dist1;
  delete Dist2;
  delete Dist3;
  delete Dist4;
  out1.close();
}


void AnalyzeStar::DoLadybot()
{

  TProfile* LadyProfile= new TProfile("LadybotProfile","Ladibot",80,0.,4000.,0.,120.,"s");
  
  
  //TCanvas *CanLadybot= new TCanvas("Ladybot","Ladybot",800,800);
  //CanLadybot->Divide(1,3);
  TCanvas *CanLadybot= new TCanvas("Ladybot","Ladybot",800,800);
  CanLadybot->SetFillColor(10);
  CanLadybot->SetFrameFillColor(0);
  CanLadybot->SetHighLightColor(2);
  CanLadybot->SetBorderSize(0);
  
  
  //vector<HitStation>::iterator sta;  
  
  cout << "toto" << endl;
  //int Size=fTree->GetEntries();
  int Size=fChain->GetEntries();

  
  for(int i=0;i<Size;i++)
    {
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      LadyProfile->Fill(fHitStation->fR_sf,(double)(fHitStation->fNbOfBinsOverThres));
    }

  
  
  LadyProfile->SetMarkerStyle(8);
  LadyProfile->SetMarkerSize(1.);
  LadyProfile->GetXaxis()->SetTickLength(0.01);
  LadyProfile->GetYaxis()->SetTickLength(0.01);
  LadyProfile->GetXaxis()->SetTitle("Distance to axis in m");
  LadyProfile->GetYaxis()->SetTitle("Nb of bins over threshold");
  
  LadyProfile->Draw();	 
  return;
    }


void AnalyzeStar::DoLadybotToT()
{

  TProfile* LadyProfile= new TProfile("LadybotProfile","Ladibot",80,0.,4000.,0.,120.,"s");
  
  
  //TCanvas *CanLadybot= new TCanvas("Ladybot","Ladybot",800,800);
  //CanLadybot->Divide(1,3);
  TCanvas *CanLadybot= new TCanvas("Ladybot","Ladybot",800,800);
  CanLadybot->SetFillColor(10);
  CanLadybot->SetFrameFillColor(0);
  CanLadybot->SetHighLightColor(2);
  CanLadybot->SetBorderSize(0);
  
  int Size=fChain->GetEntries();

  for(int i=0;i<Size;i++)
    {
      //fTree->GetEntry(i);
      fChain->GetEntry(i);
      if (fHitStation->fT2ToT)
	{
	  LadyProfile->Fill(fHitStation->fR_sf,(double)(fHitStation->fNbOfBinsOverThres));
	}
    }
  
  LadyProfile->SetMarkerStyle(8);
  LadyProfile->SetMarkerSize(1.);
  LadyProfile->GetXaxis()->SetTickLength(0.01);
  LadyProfile->GetYaxis()->SetTickLength(0.01);
  LadyProfile->GetXaxis()->SetTitle("Distance to axis in m");
  LadyProfile->GetYaxis()->SetTitle("Nb of bins over threshold");
  
  LadyProfile->Draw();	
  
  return;
}

void AnalyzeStar::DoLTP(char *energy)
{
  
  gStyle->SetOptStat(0);
  
  string tta;
  for (int j=0;j<6;j++)
    {
      

      
      if (j==0)
	tta="0.00";
      if (j==1)
	tta="25.84";
      if (j==2)
	tta="36.87";
      if (j==3)
	tta="45.57";
      if (j==4)
	tta="53.13";
      if (j==5)
	tta="60.00";
      
      TProfile* ToTProfile1 = new TProfile("T2ToTProfile1","LTP T2ToT",80,0.,4000.,0.,2.,"s");
      TProfile* T2thProfile1 = new TProfile("T2ThProfile1","LTP T2th",80,0.,4000.,0.,2.,"s");
      TProfile* T1Profile1 = new TProfile("T1Profile1","LTP T1",80,0.,4000.,0.,2.,"s");
      TProfile* T2Profile1 = new TProfile("T2Profile","LTP T2",80,0.,4000.,0.,2.,"s");
      
  
      TCanvas *CanLTP=new TCanvas("LTP","LTP",10,10,1136,772);
      
      CanLTP->SetFillColor(10);
      CanLTP->SetFrameFillColor(0);
      CanLTP->SetHighLightColor(2);
      CanLTP->SetBorderSize(0);
      CanLTP->Divide(1,4);
      CanLTP->SetFillColor(0);
      
      char * theta=new char[10];
      sprintf(theta,"%s",tta.c_str());
      ReadFilesCDF(energy,theta);
      
      cout << fChain->GetEntries() << " stations " << endl;
      for (int i=0;i< fChain->GetEntries() ; i++) 
	{

	  fChain->GetEntry(i);
      
	  if (fHitStation->fT1Threshold==1)
	    {
	      T1Profile1->Fill(fHitStation->fR_sf,1.);
	    }
	  else
	    {
	      T1Profile1->Fill(fHitStation->fR_sf,0.);
	    }
	  
	  if (fHitStation->fT2Threshold==1)
	    {
	      T2thProfile1->Fill(fHitStation->fR_sf,1.);
	    }
	  else
	    {
	      T2thProfile1->Fill(fHitStation->fR_sf,0.);
	    }
	  
	  if ((fHitStation->fT2Threshold==1) || (fHitStation->fT2ToT==1))
	    {
	      T2Profile1->Fill(fHitStation->fR_sf,1.);
	    }
	  else
	    {
	      T2Profile1->Fill(fHitStation->fR_sf,0.);
	    }
      
	  if (fHitStation->fT2ToT==1)
	    {
	      ToTProfile1->Fill(fHitStation->fR_sf,1.);
	    }
	  else
	    {
	      ToTProfile1->Fill(fHitStation->fR_sf,0.);
	    }
	
	  if (fHitStation->fT2ToT==1)
	    {
	      ToTProfile1->Fill(fHitStation->fR_sf,1.);
	    }
	  else
	    {
	      ToTProfile1->Fill(fHitStation->fR_sf,0.);
	    }  
	}
  
      
      
      CanLTP->cd(1);
      T1Profile1->SetLineColor(1);
      T1Profile1->SetMarkerColor(1);
      T1Profile1->SetMarkerStyle(8);
      T1Profile1->SetLineStyle(1);
      T1Profile1->GetXaxis()->SetTickLength(0.01);
      T1Profile1->GetYaxis()->SetTickLength(0.01);
      T1Profile1->GetXaxis()->SetTitle("Distance to axis in m");
      T1Profile1->GetYaxis()->SetTitle("Probability");
      T1Profile1->SetMinimum(0.);
      T1Profile1->SetMaximum(1.3);
      
      T1Profile1->Draw();
      
      CanLTP->cd(2);
      
      T2Profile1->SetLineColor(1);
      T2Profile1->SetMarkerColor(1);
      T2Profile1->SetMarkerStyle(8);
      T2Profile1->SetLineStyle(1);
      T2Profile1->GetXaxis()->SetTickLength(0.01);
      T2Profile1->GetYaxis()->SetTickLength(0.01);
      T2Profile1->SetMinimum(0.);
      T2Profile1->SetMaximum(1.3);
      T2Profile1->Draw();
      
      
      CanLTP->cd(3);
      
      T2thProfile1->SetLineColor(1);
      T2thProfile1->SetMarkerColor(1);
      T2thProfile1->SetMarkerStyle(8);
      T2thProfile1->SetLineStyle(1);
      T2thProfile1->GetXaxis()->SetTickLength(0.01);
      T2thProfile1->GetYaxis()->SetTickLength(0.01);
      T2thProfile1->SetMinimum(0.);
      T2thProfile1->SetMaximum(1.3);
      T2thProfile1->Draw();
  
  
  
  
      CanLTP->cd(4);
  
      ToTProfile1->SetLineColor(1);
      ToTProfile1->SetMarkerColor(1);
      ToTProfile1->SetMarkerStyle(8);
      ToTProfile1->SetLineStyle(1);
      ToTProfile1->GetXaxis()->SetTickLength(0.01);
      ToTProfile1->GetYaxis()->SetTickLength(0.01);
      ToTProfile1->SetMinimum(0.);
      ToTProfile1->SetMaximum(1.3);
      ToTProfile1->Draw();
 

      ofstream out1;
      char *writename=new char[100];
      double EgeV=atof(energy);
      int E=ceil(10*log10(EgeV*1e9));
      int TTA=(int)(atoi(tta.c_str()));
      cout << EgeV << "  " << E << endl;
      
      sprintf(writename,"ltp_p_%d_%d.txt",E,TTA);
      out1.open(writename,ios::out);
      //out1.open("ltp_p_174_0.txt",ios::out);
      double t1,t2th,t2,tot,d;
  //  double rmst1,rmst2,rmstot;
      //out1 << "Distance" << "  " << "T1" << "   " << "RMST1" << "  "  << "T2" << "  " << "RMST2" << "  "  << "ToT" << "  " << "RMSToT" << endl;
      out1 << "Distance" << "  " << "T1"  << "  "  << "T2Th" <<  "  "  << "T2ToT" << "  " << "T2" << endl;
      for (int k=1;k<=80;k++)
	{
	  t1=T1Profile1->GetBinContent(k);
	  t2=T2Profile1->GetBinContent(k);
	  t2th=T2thProfile1->GetBinContent(k);
	  tot=ToTProfile1->GetBinContent(k);
	  //rmst1=T1Profile1->GetBinError(k);
	  //rmst2=T2Profile1->GetBinError(k);
	  //rmstot=ToTProfile1->GetBinError(k);
	  d=50.*(k-1.);
	  //cout << d << "  " << t1 << "  " <<  rmst1 << "  " <<t2 << "  " << rmst2 << "  " <<tot  << "  " << rmstot << endl;
	  //out1 << d << "  " << t1 << "  " <<  rmst1 << "  " <<t2 << "  " << rmst2 << "  " <<tot << "  " << rmstot << endl;  
	  out1 << d << "  " << t1 << "  "  << t2th << "  " << tot << "  " << t2 << endl;
	}
      out1.close();
      
    }
 return;
}

//    out1.open("ltp_182_45.txt",ios::out);
//    out1 << "Distance" << "  " << "T1" << "   " << "RMST1" << "  "  << "T2" << "  " << "RMST2" << "  "  << "ToT" << "  " << "RMSToT" << endl;
//    for (int k=1;k<=80;k++)
//      {
//        t1=T1Profile2->GetBinContent(k);
//        t2=T2Profile2->GetBinContent(k);
//        tot=ToTProfile2->GetBinContent(k); 
//        rmst1=T1Profile2->GetBinError(k);
//        rmst2=T2Profile2->GetBinError(k);
//        rmstot=ToTProfile2->GetBinError(k);
//        d=50.*(k-1.);cout << d << "  " << t1 << "  " <<  rmst1 << "  " <<t2 << "  " << rmst2 << "  " <<tot  << "  " << rmstot << endl;
//        out1 << d << "  " << t1 << "  " <<  rmst1 << "  " <<t2 << "  " << rmst2 << "  " <<tot << "  " << rmstot << endl;

   
      
//      }
//    out1.close();
  
//    out1.open("ltp_182_60.txt",ios::out);
//    out1 << "Distance" << "  " << "T1" << "   " << "RMST1" << "  "  << "T2" << "  " << "RMST2" << "  "  << "ToT" << "  " << "RMSToT" << endl;
//    for (int k=1;k<=80;k++)
//      {
//        t1=T1Profile3->GetBinContent(k);
//        t2=T2Profile3->GetBinContent(k);
//        tot=ToTProfile3->GetBinContent(k);
//        rmst1=T1Profile3->GetBinError(k);
//        rmst2=T2Profile3->GetBinError(k);
//        rmstot=ToTProfile3->GetBinError(k);
//        d=50.*(k-1.);
//        cout << d << "  " << t1 << "  " <<  rmst1 << "  " <<t2 << "  " << rmst2 << "  " <<tot  << "  " << rmstot << endl;
//        out1 << d << "  " << t1 << "  " <<  rmst1 << "  " <<t2 << "  " << rmst2 << "  " <<tot << "  " << rmstot << endl;

      
      
//      }
//    out1.close();



void AnalyzeStar::DrawSaturatedTrace()
{
  //short adcsum[2][MAXNUMBEROFADCBINS];
  //short adcsum_mu[2][MAXNUMBEROFADCBINS];
  //short adcsum_ns[MAXNUMBEROFADCBINS];
  TH1F* histo;
  TH1F* histo_ns;
  TH1F* histo_mu;
  
  int nmax=MAXNUMBEROFADCBINS;
  
  /*for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    {
      adcsum_ns[i]=0;
      for(int j = 0;j<2;j++)
	{
	  adcsum[j][i]=0;
	  adcsum_mu[j][i]=0;
	}
	}*/
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  fCanvasList->Add(new TCanvas("sat","saturated trace ",400,300));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  int idistance = (int)fHitStation->fR_sf;
  int id = fHitStation->fId;
  
  
  char*hname = new char[15];
  sprintf(hname,"SUMADC_%04d_%04d",id,idistance);
  
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,400,0.,400);
  delete hname;
  
  hname = new char[18]; 
  sprintf(hname,"SUMADC_%04d_%04d_ns",id,idistance);
  cout<<hname<<endl;
  histo_ns=  new TH1F(hname,hname,400,0.,400);
  delete hname;
  
  hname = new char[18];
  sprintf(hname,"SUMADC_%04d_%04d_mu",id,idistance);
  
  cout<<hname<<endl;
  histo_mu=  new TH1F(hname,hname,400,0.,400);
  delete hname;
  
  
 // adcsum=fHitStation->fADC;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fADC[1][k]!=0){
      histo->Fill(k,fHitStation->fADC[1][k]);
    }
  
  //adcsum_mu=fHitStation->fADC_mu;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fADC_mu[1][k]!=0){
      histo_mu->Fill(k,fHitStation->fADC_mu[1][k]);
    }
  
 
  // adcsum_ns=fHitStation->fADC_ns;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fADC_ns[k]!=0){
      histo_ns->Fill(k,fHitStation->fADC_ns[k]);
    }
  
  
  histo->SetLineColor(kRed);
  histo_ns->SetLineColor(kCyan);
  histo_mu->SetLineColor(kBlue);
  histo_ns->Draw();
  histo->Draw("SAME");
  histo_mu->Draw("SAME");
  
  
  
}

void AnalyzeStar::DrawSaturatedTrace(int ipm)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  //short adcsum_ns[MAXNUMBEROFADCBINS];
  //short adcsum_mu[2][MAXNUMBEROFADCBINS];
  TH1F* histo;
  TH1F* histo_ns;
  TH1F* histo_mu;
 
  
  
  int nmax=MAXNUMBEROFADCBINS;
  
  /* for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
    {
       adcsum_ns[i]=0;
       for(int j = 0;j<2;j++)
	 {
	   adcsum[j][i]=0;
	   adcsum_mu[j][i]=0;
	 }
	 }*/
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  
  fCanvasList->Add(new TCanvas("sat","saturated trace ",400,300));
  fCanvas=(TCanvas *) fCanvasList->Last();
  
  
  int idistance = (int)fHitStation->fR_sf;
  int id = fHitStation->fId;
  
  char*hname ;
  hname = new char[17];
  sprintf(hname,"ADC_PM%1d_%04d_%04d",ipm,id,idistance);
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,400,0.,400);
  delete hname;
  
  hname = new char[21]; 
  sprintf(hname,"ADC_PM%1d_%04d_%04d_mu",ipm,id,idistance);
  cout<<hname<<endl;
  histo_mu=  new TH1F(hname,hname,400,0.,400);
  cout<<"histo creee"<<endl;
  delete hname;
  
  
  hname = new char[21]; 
  sprintf(hname,"ADC_PM%1d_%04d_%04d_ns",ipm,id,idistance);
  cout<<hname<<endl;
  histo_ns=  new TH1F(hname,hname,400,0.,400);
  cout<<"histo creee"<<endl;
  delete hname;
  
  // adcsum=fHitStation->fPMT[ipm-1].fADC;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fPMT[ipm-1].fADC[1][k]!=0){
      histo->Fill(k,fHitStation->fPMT[ipm-1].fADC[1][k]);
    }
  
  //adcsum_mu=fHitStation->fPMT[ipm-1].fADC_mu;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fPMT[ipm-1].fADC_mu[1][k]!=0){
      histo_mu->Fill(k,fHitStation->fPMT[ipm-1].fADC_mu[1][k]);
    }
  
  
  //adcsum_ns=fHitStation->fPMT[ipm-1].fADC_ns;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fPMT[ipm-1].fADC_ns[k]!=0){
      histo_ns->Fill(k,fHitStation->fPMT[ipm-1].fADC_ns[k]);
    }
  
  
  histo->SetLineColor(kRed);
  histo_ns->SetLineColor(kCyan);
  histo_mu->SetLineColor(kBlue);
  histo_ns->Draw();
  histo->Draw("SAME");
  histo_mu->Draw("SAME");
  
  
}



void AnalyzeStar::GetStation(Int_t sta)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  int nmax=MAXNUMBEROFADCBINS;
 
  if (sta>= fTree->GetEntries()) 
    {
      cout<<"Error in AnalyzeStar::GetShower() ; only " << (Int_t) fTree->GetEntries() << " events in Tree"<<endl;
      return;  
    }	
 
 
   cout<<"test"<<endl;
   fTree->GetEvent(sta);
  
 
  
   cout << "\t tank " <<fHitStation->fId <<"\t pel = "<<fHitStation->fNpe<<
     "\t nmu = "<<fHitStation->fNmu<<"\t dist = "<<fHitStation->fR_sf<<endl ; 
 
   if(fHitStation->fSampFact!=1.)cout<<"careful sampling factor applied to save CPU = " << fHitStation->fSampFact<<endl;
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  TH1F* histo;
  char*hname = new char[15];
  int idistance = (int)fHitStation->fR_sf;
  int id = fHitStation->fId;
        
  sprintf(hname,"SUMADC_%04d_%04d",id,idistance);
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,400,0.,400);
  delete hname;   
  //adcsum=fHitStation->fADC;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fADC[0][k]!=0){
      histo->Fill(k,fHitStation->fADC[0][k]);
    }


  
  histo->SetLineColor(kRed);
  histo->Draw();
} 

void AnalyzeStar::GetStation(Int_t sta,Int_t gain)
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  int nmax=MAXNUMBEROFADCBINS;
 
  if (sta>= fTree->GetEntries()) 
    {
      cout<<"Error in AnalyzeStar::GetShower() ; only " << (Int_t) fTree->GetEntries() << " events in Tree"<<endl;
      return;  
    }	
 
 
   cout<<"test"<<endl;
   fTree->GetEvent(sta);
  
 
  
   cout << "\t tank " <<fHitStation->fId <<"\t pel = "<<fHitStation->fNpe<<
     "\t nmu = "<<fHitStation->fNmu<<"\t dist = "<<fHitStation->fR_sf<<endl ; 
 
   if(fHitStation->fSampFact!=1.)cout<<"careful sampling factor applied to save CPU = " << fHitStation->fSampFact<<endl;
  
  gStyle->SetOptStat(1000000); 
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleH(0.08);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  TH1F* histo;
  char*hname = new char[15];
  int idistance = (int)fHitStation->fR_sf;
  int id = fHitStation->fId;
        
  sprintf(hname,"SUMADC_%04d_%04d",id,idistance);
  cout<<hname<<endl;
  histo=  new TH1F(hname,hname,400,0.,400);
  delete hname;   
  //adcsum=fHitStation->fADC;
  for(int k=0; k<nmax; k++)
    if(fHitStation->fADC[gain][k]!=0){
      histo->Fill(k,fHitStation->fADC[gain][k]);
    }


  
  histo->SetLineColor(kRed);
  histo->Draw();
} 



void AnalyzeStar::LateralProfile()
{
  // short adcsum[2][MAXNUMBEROFADCBINS];
  //short adcsum_ns[MAXNUMBEROFADCBINS];
  int nmax = MAXNUMBEROFADCBINS;
  int sumtot,sumsat;

  int nsta=(int)fTree->GetEntries();
  TProfile* latsat = new TProfile("latsat","latsat",2500,0,2500,0,500000); 
  TProfile* lat = new TProfile("lat","lat",2500,0,2500,0,500000); 
  
 
  
  for(int sta=0;sta<nsta;sta++)
    {
      
      sumtot=0;
      sumsat=0;
      //for(int i = 0;i<MAXNUMBEROFADCBINS;i++)
      //{
      //  adcsum_ns[i]=0;
      //  for(int j = 0;j<2;j++)
      //  {
      //  adcsum[j][i]=0;
	  
      //  }
      //  }
      
      
      fTree->GetEvent(sta);
      
      
      
      cout << "\t tank " <<fHitStation->fId <<"\t pel = "<<fHitStation->fNpe<<
	"\t nmu = "<<fHitStation->fNmu<<"\t dist = "<<fHitStation->fR_sf<<endl ; 
      
      if(fHitStation->fSampFact!=1.)cout<<"careful sampling factor applied to save CPU = " << fHitStation->fSampFact<<endl;
      
      
      // adcsum=fHitStation->fADC;
      for(int k=0; k<nmax; k++)
	sumsat+=fHitStation->fADC[1][k];
      latsat->Fill(fHitStation->fR_sf,sumsat);
      
      //adcsum_ns=fHitStation->fADC_ns;
      for(int k=0; k<nmax; k++)
	sumtot+=fHitStation->fADC_ns[k];
      lat->Fill(fHitStation->fR_sf,sumtot);
      cout<<sumsat<<" "<<sumtot<<endl;
    }

 lat->SetMarkerColor(kRed);
 latsat->SetMarkerColor(kCyan);
  
 lat->Draw();
 latsat->Draw("SAME");
  
}



void AnalyzeStar::PlotNbotDistr()
{

  double dist1=300.;
  double dist2=600.;
  double dist3=1000.;
  double dist4=1500.;
  double delta=100.;
  TCanvas *Nbots= new TCanvas("Distribution","Distribution of bins over threshold",10,10,1136,772);

  Nbots->SetFillColor(10);
  Nbots->SetFrameFillColor(0);
  Nbots->SetHighLightColor(2);
  Nbots->SetBorderSize(0);
  Nbots->Divide(2,2);

  
  TH1D *Dist1 = new TH1D("Dist1","d=300m",100,0.,60);
  TH1D *Dist2 = new TH1D("Dist2","d=600m",100,0.,60);
  TH1D *Dist3 = new TH1D("Dist3","d=1000m",100,0.,60);
  TH1D *Dist4 = new TH1D("Dist4","d=1500m",100,0.,60);
   
  Dist1->GetXaxis()->SetTitle("Nb of bins over threshold");
  Dist1->GetYaxis()->SetTitle("Count");
  Dist1->GetXaxis()->SetLabelSize(0.05);
  Dist1->GetXaxis()->SetTitleSize(0.045);
  Dist1->GetXaxis()->SetTitleOffset(1.1);
  Dist1->GetYaxis()->SetLabelSize(0.05);
  Dist1->GetYaxis()->SetTitleSize(0.045);
  Dist1->GetYaxis()->SetTitleOffset(1.);

  Dist2->GetXaxis()->SetTitle("Nb of bins over threshold");
  Dist2->GetYaxis()->SetTitle("Count");
  Dist2->GetXaxis()->SetLabelSize(0.05);
  Dist2->GetXaxis()->SetTitleSize(0.045);
  Dist2->GetXaxis()->SetTitleOffset(1.1);
  Dist2->GetYaxis()->SetLabelSize(0.05);
  Dist2->GetYaxis()->SetTitleSize(0.045);
  Dist2->GetYaxis()->SetTitleOffset(1.);


  Dist3->GetXaxis()->SetTitle("Nb of bins over threshold");
  Dist3->GetYaxis()->SetTitle("Count");
  Dist3->GetXaxis()->SetLabelSize(0.05);
  Dist3->GetXaxis()->SetTitleSize(0.045);
  Dist3->GetXaxis()->SetTitleOffset(1.1);
  Dist3->GetYaxis()->SetLabelSize(0.05);
  Dist3->GetYaxis()->SetTitleSize(0.045);
  Dist3->GetYaxis()->SetTitleOffset(1.);

  Dist4->GetXaxis()->SetTitle("Nb of bins over threshold");
  Dist4->GetYaxis()->SetTitle("Count");
  Dist4->GetXaxis()->SetLabelSize(0.05);
  Dist4->GetXaxis()->SetTitleSize(0.045);
  Dist4->GetXaxis()->SetTitleOffset(1.1);
  Dist4->GetYaxis()->SetLabelSize(0.05);
  Dist4->GetYaxis()->SetTitleSize(0.045);
  Dist4->GetYaxis()->SetTitleOffset(1.);





  
  for (int i=0;i< fChain->GetEntries() ; i++) 
    {
      fChain->GetEntry(i);
      //fStList = fEvent->fHitStationList;
      
      //for(sta=fStList.begin();  sta!= fStList.end(); sta++)
      //{ 
      if (fHitStation->fT2ToT)
	{
	  
	  if ( (fHitStation->fR_sf > (dist1-delta)) && (fHitStation->fR_sf< (dist1+delta)) )
	    {
	      Dist1->Fill((double)(fHitStation->fNbOfBinsOverThres));
	    }
	  
	  if ((fHitStation->fR_sf> (dist2-delta)) && (fHitStation->fR_sf< (dist2+delta)))
	    {
	      Dist2->Fill((double)(fHitStation->fNbOfBinsOverThres));
	    }
	  
	  if ((fHitStation->fR_sf> (dist3-delta)) && (fHitStation->fR_sf< (dist3+delta)))
	    {
	      Dist3->Fill((double)(fHitStation->fNbOfBinsOverThres));
	    } 
	  
	  if ((fHitStation->fR_sf> (dist4-delta)) && (fHitStation->fR_sf< (dist4+delta)))
	    {
	      Dist4->Fill((double)(fHitStation->fNbOfBinsOverThres));
	    } 
	}
      
    }     
  

  Nbots->cd(1);
  Dist1->Draw();
  Nbots->cd(2);
  Dist2->Draw();	
  Nbots->cd(3);
  Dist3->Draw();	
  Nbots->cd(4);
  Dist4->Draw();	
  
  double mean1=Dist1->GetMean();
  double rms1=Dist1->GetRMS();
  double mean2=Dist2->GetMean();
  double rms2=Dist2->GetRMS();
  double mean3=Dist3->GetMean();
  double rms3=Dist3->GetRMS();
  double mean4=Dist4->GetMean();
  double rms4=Dist4->GetRMS();
  cout << "mean1= " << mean1 << "  rms1="  << rms1 << endl;
  cout << "mean2= " << mean2 << "  rms2="  << rms2 << endl; 
  cout << "mean3= " << mean3 << "  rms3="  << rms3 << endl;
  cout << "mean4= " << mean4 << "  rms4="  << rms4 << endl;
  
}

 void AnalyzeStar::PrintS1000(char *energy)
{
 ofstream out1;  
 out1.open("energyS1000.out",ios::app);
 gStyle->SetOptStat(0);
 
 string tta;
  for (int j=0;j<6;j++)
    {
      if (j==0)
	tta="0.00";
      if (j==1)
	tta="25.84";
      if (j==2)
	tta="36.87";
      if (j==3)
	tta="45.57";
      if (j==4)
	tta="53.13";
      if (j==5)
	tta="60.00";
     
      
      char * theta=new char[10];
      sprintf(theta,"%s",tta.c_str());
      ReadFilesCDF(energy,theta);
      double s1000=0;
      int ns1000=0;
      cout << fChain->GetEntries() << " stations " << endl;
      double EgeV=atof(energy);
      double E=ceil(10*log10(EgeV*1e9));
      int TTA=(int)(atoi(tta.c_str()));
      cout << EgeV <<"  " <<E<<" "<<TTA<< endl;
	  
      for (int i=0;i< fChain->GetEntries() ; i++) 
	{
	  
	  fChain->GetEntry(i);
      
	  
	  if(fHitStation->fR_sf>980 && fHitStation->fR_sf<1020){
	    s1000+=fHitStation->fSignalInVem;
	    ns1000++;
	    cout<<fHitStation->fR_sf<<" "<<ns1000<<endl;
	  }
	
	}
    
      
      out1 << E << "  " <<TTA<<" "<<s1000/(double)ns1000 <<endl;
	  
	
	 
    }
  
  
  out1.close();
  
  return;
  


  
}



void AnalyzeStar::PrintTree()
{  
  fTree->Print();
}

void AnalyzeStar::ReadFile(char * filename,int opt)
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
  
  
  
  
  
  fBranch= fTree->GetBranch("HitStation");
  fTree->SetBranchAddress("hstation",&fHitStation);
  
  cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
  
  fTree->GetEntry(0);
  
  return;
}
//---------------------------------------------------------------------



void AnalyzeStar::ReadFilesCDF(char* Energy, char * tta)
{
  
  if(fChain!=0)
    { 
      delete fChain;
      fChain=0;
    }
  fChain= new TChain("TreeEvent");
  
  ifstream in1;
  // string Sim="/projet/auger/Corsika/CDF/";
  string  Primary="proton/";
  //string  Theta="0.00d/";
  // 0.00  25.84 36.87  45.57  53.13  60.00

  //  1.585e+09/ 1.585e+11/ 1e+10/  2.512e+08/ 2.512e+10/ 3.981e+08/ 3.981e+10/ 6.31e+09/
//  1.585e+10/ 1e+09/  1e+11/ 2.512e+09/ 2.512e+11/ 3.981e+09/ 6.31e+08/  6.31e+10/
  //string Energy="2.512e+08/";
  string DebutFichier="Sim_stararray_DAT00*";
  string tampon="tampon.dat";  
  int Phi;
  string commandrm="rm tampon.dat";
  char * readfilename1= new char[250];
  char * readfilename2= new char[250];; 
  char * commandls= new char[250];;

  for (int i=0;i<5;i++)
    {
     
      if (i==0)
	Phi=54;
      if (i==1)
	Phi=126;
      if (i==2)
	Phi=198;
      if (i==3)
	Phi=270;
      if (i==4)
	Phi=342;

      //sprintf(commandls,"ls %s%s%sd/%s%dd/%s > %s",Sim.c_str(),Primary.c_str(),Theta.c_str(),Energy.c_str(),Phi,DebutFichier.c_str(),tampon.c_str());
      
      sprintf(commandls,"ls %s%s%sd/%s/%dd/%s > %s",Sim.c_str(),Primary.c_str(),tta,Energy,Phi,DebutFichier.c_str(),tampon.c_str());
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
  
  fBranch= fChain->GetBranch("HitStation");
  fChain->SetBranchAddress("hstation",&fHitStation);
  
  cout << " \t * There is " << fChain->GetEntries() << " events in the input file "<<endl;
  
  fChain->GetEntry(0);
  fChain->Print();
  return;

//    //fFile = new TFile(filename,"READ");	// Open file	
//    //fFile->ls();
  
//    //fTree = (TTree *) fFile->Get("TreeEvent");      // Get Tree


//    fChain->Print();   
//    //fTree = (TTree *) fTree->Get("TreeEvent");   
//    //fTree->Print(); 
  
//    if(fBranch!=0)
//      { 
//        delete fBranch;
//        fBranch=0;
//      }
  
  
//    fChain->SetBranchAddress("event",&fEvent);
  
  
//    //fBranch= fTree->GetBranch("event");
//    //fTree->SetBranchAddress("event",&fEvent);
  
//    //cout << " \t * There is " << fTree->GetEntries() << " events in the input file "<<endl;
//    cout << " \t * There is " << fChain->GetEntries() << " events in the input file "<<endl;
//    fChain->ls();
//    //fTree->GetEntry(0);
//    fChain->GetEntry(0);
  
//    return;
}


//---------------------------------------------------------------------
void AnalyzeStar::ReadFilesLTP()
{
  cout << "toto" << endl;  
  if(fChain!=0)
    { 
      delete fChain;
      fChain=0;
    }
  fChain= new TChain("TreeEvent");
  cout << "toto" << endl;  
 
 
  string file1("/projet/auger/Corsika/CDF/proton/0.00d/1e+09/126d/Sim_stararray_DAT001183.root");
  string file2("/projet/auger/Corsika/CDF/proton/0.00d/1e+09/126d/Sim_stararray_DAT000075.root");
  
  
  
  //sprintf(readfilename,"%s_%s_%d_%d_%d.root",Sim.c_str(),Primary.c_str(),energy,tta,i); 
  cout << "ajout du fichier " << file1 << endl;
  fChain->Add(file1.c_str());
  cout << "ajout du fichier " << file2 << endl;
  fChain->Add(file2.c_str());
  
  if(fBranch!=0)
    { 
      delete fBranch;
      fBranch=0;
    }
   
  fBranch= fChain->GetBranch("HitStation");
  fChain->SetBranchAddress("hstation",&fHitStation);
  
  cout << " \t * There is " << fChain->GetEntries() << " events in the input file "<<endl;
  
  fChain->GetEntry(0);
  fChain->Print();
  return;
}







  




