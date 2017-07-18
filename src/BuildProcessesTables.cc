
/*
  BuildProcessesTables.cc

  implementation file for 

*/

using namespace std;

#include "BuildProcessesTables.h"
#include "Calib.h"
#include "Constants.h"
#include "EasySim.h"
#include <math.h>
#include <iostream>


double gCherFact;
double gCherPrb[NWAVEL];
double gTMinForDeltaRays;
double gBremTable[NSTEP_E][NDIV];
double gRepartFuncMichelSpectrum[NENER_ELEC_MUDECAY];
double gMuonEnergyAngularSpectrum[ANGULAR_FLUX_SIZE];
double gMuonProbaAngularSpectrum[ANGULAR_FLUX_SIZE];

double Integrande1(double ener)
{
  double fact = pow(GFERMI,2.)*pow(MASS_MU, 5.)/192./pow(PI, 3.);
  double e = ener/(MASS_MU/2.);
  
  return fact*2*e*e*(3.-2*e);
}

void BuildProcessesTables()
{
  double term,cherprb[NWAVEL];
  for (int i=0;i<NWAVEL;i++)
    {
      term= PMTQE[i]/pow((WAVELMIN+(i+.5)*DWAVEL),2);
      // term= 1/pow((WAVELMIN+(i+.5)*DWAVEL),2);
      if(i==0) cherprb[i]=term;
      else cherprb[i]=cherprb[i-1]+term;         
    }
  
  gCherFact=CHERCONST*cherprb[NWAVEL-1]*DWAVEL*COLLEF;
  //gCherFact=CHERCONST*cherprb[NWAVEL-1]*DWAVEL;
  
  for(int i=0;i<NWAVEL;i++){
    gCherPrb[i] =cherprb[i]/cherprb[NWAVEL-1];
    //  cout<<"88888 "<<(WAVELMIN+(i+.5)*DWAVEL)*1.e9<<" "<<gCherPrb[i]<<endl;
  }
  
  
  //calculation of TMIN for delta rays
  gTMinForDeltaRays = MASS_E*(1/sqrt(1-1./(WATINDEX*WATINDEX))-1);

  //*****************************************************************
  //make Bremstrahlung tables 
  
  /* This and following function make a look-up table for bremsstrahlung
     effect in electrons. The differential prob. per g cm**-2 of an interaction
     which transfers 1 to 99% of the incident electron K.E is calculated using  
     functions given in Rossi page 49. The resulting surface is then
     integrated and normalised to make an efficient table look-up easy. */ 
    
  int     ebin,i;
  double  T,U,Nu,MaxNu,dNu,PhiOxygen,PhiHydrogen,Phi;
  double  IntPhi;
  
  
  cerr<<"\tBuilding Bremsstrahlung look-up table"<<endl;
  
  
  for(ebin=0;ebin<NSTEP_E;ebin++)
    {
      T=ESTEP_E[ebin];
      //T=0.1*pow(10,(ebin+0.5)*0.2);
      U=T+MASS_E;
      MaxNu=1-(MASS_E/U);
      dNu=MaxNu/NDIV;
      IntPhi=0;
      
      for(i=1;i<NDIV;i++)
	{ 
	  Nu=(i)*dNu;
	  PhiOxygen  =CalcPhiUNudNu(U,Nu,dNu,15.999,8.0);
	  PhiHydrogen=CalcPhiUNudNu(U,Nu,dNu, 1.008,1.0);
	  Phi=PhiOxygen*(16.0/18)+PhiHydrogen*(2.0/18);
	  IntPhi+=Phi;
	  gBremTable[ebin][i]=IntPhi;
	}
      
      gBremTable[ebin][0]=IntPhi;
    }
  
  /* Normalise the table */
  for(ebin=0;ebin<NSTEP_E;ebin++)
    {
      for(i=1;i<NDIV;i++)
	{
	  gBremTable[ebin][i]/=gBremTable[ebin][0];
	  //cout<<i<<" "<<ebin<<" "<<gBremTable[ebin][i]<<endl;
	} 
    }
  
  
  //*****************************************************************
  // make table for energy electron distribution after muon decay
  for (int i=0; i<NENER_ELEC_MUDECAY; i++)
    double ener = MASS_E*pow(10, i*log10(MASS_MU/2./MASS_E)/(double)(NENER_ELEC_MUDECAY-1.));      
  //      gRepartFuncMichelSpectrum[i] = qromb(&Integrande1, 0., ener);
  //cout << ener << "  " << Integrande1(ener) << endl;
  
  
  
  //build muons spectrum according to zenital angle 
  if((theConfig()->PartEnergy==0)&&(theConfig()->PartCode==3) && 
     ((theConfig()->PartMode=="SCINTILLATOR")||(theConfig()->PartMode=="FIXEDTHETA")||(theConfig()->PartMode=="VEM")||(theConfig()->PartMode=="HORIZONTAL"))){
   cerr<<"\tBuilding muon energy integrated distribution as a function of angle look-up table"<<endl;
   int NBins=60000;
   int NBinsPerDecade=10000;
   int NNewBins= ANGULAR_FLUX_SIZE;
   int NNewBinsPerDecade= NNewBins*NBinsPerDecade/NBins ;
   
   double zenith,maxval;
   double C,g0,g1,g2,g3,g4;
   double prob[ANGULAR_FLUX_SIZE];
   double y[ANGULAR_FLUX_SIZE];
   if(theConfig()->PartMode=="FIXEDTHETA") 
     zenith =theConfig()->PartTheta;
   else if(theConfig()->PartMode=="SCINTILLATOR") {
     Double_t xhigh= XPOS_SCINT_HIGH;
     Double_t yhigh= YPOS_SCINT_HIGH;
     Double_t xlow= XPOS_SCINT_LOW;
     Double_t ylow= YPOS_SCINT_LOW;
     Double_t rad = sqrt((xhigh-xlow)*(xhigh-xlow)+(yhigh-ylow)*(yhigh-ylow));
     if(ZPOS_SCINT_HIGH!=ZPOS_SCINT_LOW)
	zenith = atan(rad/(ZPOS_SCINT_HIGH-ZPOS_SCINT_LOW));
      else 
	zenith=PI/2.-0.1;
    }
    else if (theConfig()->PartMode=="VEM")
      {
	zenith=0.;
      }
    else if (theConfig()->PartMode=="HORIZONTAL")
      {
	zenith=PI/2-0.1;
      }
    else {
      cerr<<" Error in definition of mode"<<endl; 
    }
    
    for(int i=0;i<NBins;i++){
      double en=pow(10,(double)(i)/(double)NBinsPerDecade);
      double P=sqrt(en*en-0.1057*0.1057);
      if( P <= 9.2765e2 ) {
	
	C = 2.950e-3;
	g0 = 0.3061;
	g1 = 1.2743;
	g2 = -0.2630;
	g3 = 0.0252;
	
      } else if( P > 9.2765e2 && P <= 1.5878e3 ) {
	
	C = 1.781e-2;
	g0 = 1.7910;
	g1 = 0.3040;
	g2 = 0;
	g3 = 0;
	
      } else if( P > 1.5878e3 && P <= 4.1625e5 ) {
	
	C = 1.435e1;
	g0 = 3.6720;
	g1 = 0;
	g2 = 0;
	g3 = 0;
	
      } else {
	
	C = 1.e3;
	g0 = 4;
	g1 = 0;
	g2 = 0;
	g3 = 0;
	
      }
      
      double expo = g0 + g1*log10(P) +
	g2*log10(P)*log10(P) + g3*log10(P)*log10(P)*log10(P); 
      double Dth = C*P*pow(P, -expo);
      
      
      //Flujo experimental co angulo a 0deg.
      double cth = 1;
      
      double a = 451;
      double b = 77.2;
      double c = 9.2;
      double d = 19.8;
      
      double f1 = a/(P*cth + b)*pow(5*P + c/cth,-2.57);
      double f2 = (P + d)/(P + d/cth);
      
      double Dexp0 = P*f1*f2;
      
      
      //Flujo experimental co angulo deg.
      cth = cos(zenith);
      
      double f11 = a/(P*cth + b)*pow(5*P + c/cth,-2.57);
      double f22 = (P + d)/(P + d/cth);
      
      double Dexp_cth = P*f11*f22;
      double proba =  (Dexp_cth/Dexp0)*Dth;
      int binratio=NBinsPerDecade / NNewBinsPerDecade ;
      int j=i/binratio;
      //  cout<<"binratio "<<binratio<<" "<<i<<" "<<j<<endl;
      //  cout<<i<<" "<<endl;     
      if(i-j*binratio==0){
	
	gMuonEnergyAngularSpectrum[j]=en;
	if(j==0)y[i]=0;
	else y[j]=y[j-1]+proba*
	       (gMuonEnergyAngularSpectrum[j]-gMuonEnergyAngularSpectrum[j-1]);
	if(j==(NNewBins-1))maxval=y[j];
	//	cout<<i<<" "<<j<<" "<<en<<" "<<proba<<endl;     
     }
      
    } 
    
    for(int j=1;j<NNewBins;j++){
      gMuonProbaAngularSpectrum[j]=y[j]/maxval;
     
    }
    
  }
}  
     
     
   
     
double  CalcPhiUNudNu(double U, double Nu, double dNu, double A, double Z)
{
#define ALPHA (1/137.0)
#define NAVAGADRO 6.02e+23
#define R0 2.82e-13

  double  FuncConst,GammaScreen,PhiUNudNu=0.;
  double  term1,term2,f1gam,f2gam,cgam;

  FuncConst=4*ALPHA*(NAVAGADRO/A)*(pow(Z,2)+Z)*pow(R0,2);     
  GammaScreen=100*(MASS_E/U)*(Nu/(1-Nu))*pow(Z,-0.33);
  
  if(GammaScreen<0.001)
  {
    term1=(1+pow(1-Nu,2)-(2./3)*(1-Nu))*log(183*pow(Z,-(1./3)));
    term2=(1./9)*(1-Nu);
    PhiUNudNu=FuncConst*(dNu/Nu)*(term1+term2);
  }
  if(GammaScreen>=0.001&&GammaScreen<2)
  {
    if(GammaScreen<0.8)
    {
      f1gam=-3.52*GammaScreen+20.79;
      f2gam=-2.86*GammaScreen+20.29;
    }
    else
    {
      f1gam=f2gam=-1.91*GammaScreen+19.40;
    }
    
    term1= (1+pow(1-Nu,2))*((f1gam/4)-(1./3)*log(Z));
    term2=-((2./3)*(1-Nu))*((f2gam/4)-(1./3)*log(Z));
    PhiUNudNu=FuncConst*(dNu/Nu)*(term1+term2);
  }
  if(GammaScreen>=2&&GammaScreen<=15)
  {
    cgam=0.5*exp(-0.45*GammaScreen)+0.01;
    term1=1+pow(1-Nu,2)-(2./3)*(1-Nu);
    term2=log((2*U/MASS_E)*((1-Nu)/Nu))-(1./2)-pow(cgam,1);
    PhiUNudNu=FuncConst*(dNu/Nu)*term1*term2;
  }
  if(GammaScreen>15)
  {
    term1=1+pow(1-Nu,2)-(2./3)*(1-Nu);        
    term2=log((2*U/MASS_E)*((1-Nu)/Nu))-(1./2);
    PhiUNudNu=FuncConst*(dNu/Nu)*term1*term2;
  }
  if(GammaScreen<0)  PhiUNudNu=0;  
  
  return(PhiUNudNu);
}
