/*
  Particle.cc
  implementation file for all classes Particle
*/


#include "Particle.h"
#include "EasySimConfig.h"
#include "BuildProcessesTables.h"
#include "Utils.h"
#include "Constants.h"
#include <math.h>
#include <stdlib.h>

ClassImp(Particle)
ClassImp(ParticleInTank)

extern double gCherFact;
extern double gCherPrb[NWAVEL];
extern double gNpesamples;
extern Int_t gNshapezones;
extern Int_t gNshape;
extern Int_t gShapezone[MAXNUMBEROFADCBINS];
extern double gIntshape[MAXNUMBEROFADCBINS];
extern Int_t gNbins;
extern double gTMinForDeltaRays;
extern double gBremTable[NSTEP_E][NDIV];

enum
  {
    STOP = 0,
    CONTINUE = 1,
    LAST
  };

//----------------------------------------------------------------
/*
  class Particle (parent class)
  A particle at the ground level.
  This particle has be produced by an Air Shower simulation code.
*/
//----------------------------------------------------------------


Particle::Particle()
{

}

Particle::Particle(Int_t idpart,Int_t ntop,Int_t nside,Float_t time_plane,
		   Float_t timepart,Float_t energy,Float_t cx,Float_t cy, Float_t cz )
{
  fId = idpart;
  fUX= cx;
  fUY= cy;
  fUZ= cz;
  fT= timepart;
  fT_plane= time_plane;
  fE=energy;
  fWtop = ntop;
  fWside = nside;
}

Particle::~Particle()
{

}


//----------------------------------------------------------------
/*
  class ParticleInTank
  Last modification: July 2nd, 2004
  For Detailed Mode:
  ------------------
  * Muon Simulation:
  Cherenkov production
  Delta Rays
  Muon Decay

  * Electron Simulation:
  Cherenkov production
  MultipleScattering
  Bremsstrahglung
  Ionization lost

  * Gamma Simulation:
  Compton Scattering
  Pair Production

  For Fast Mode:
  --------------
  Only Cherenkov production but NO WORKING !!!!!

*/
//----------------------------------------------------------------


ParticleInTank::ParticleInTank()
{

}

ParticleInTank::ParticleInTank (Int_t ist,Int_t id,double x,double y,double z, 
				double cx,double cy,double cz,double tim,double en,double weight)
  //  : Particle(id, 0,0,0, tim, en, cx, cy, cz)
{
  fId = id;
  fUX= cx;
  fUY= cy;
  fUZ= cz;
  fT= tim;
  //  fT_plane= time_plane;
  fE=en;

  fISta = ist;
  fX = x;
  fY = y;
  fZ = z;
  fW = weight;
  fNCherPhot = 0;
  fTCherPhot.clear();
  fNReflexions.clear();
  fIpm.clear();
  fBinContent.clear();

}
ParticleInTank::ParticleInTank (Int_t ist,Int_t id,Int_t motherid,double x,double y,double z, 
				double cx,double cy,double cz,double tim,double en,double weight)
  //  : Particle(id, 0,0,0, tim, en, cx, cy, cz)
{
  fId = id;
  fMotherId =  motherid;
  fUX= cx;
  fUY= cy;
  fUZ= cz;
  fT= tim;
  //  fT_plane= time_plane;
  fE=en;

  fISta = ist;
  fX = x;
  fY = y;
  fZ = z;
  fW = weight;
  fNCherPhot = 0;
  fTCherPhot.clear();
  fNReflexions.clear();
  fIpm.clear();
  fBinContent.clear();

}

ParticleInTank::~ParticleInTank()
{

}


// ------------------------------------------------------
// Detector simulation for each Particle
// ------------------------------------------------------

// ---------------------------
/*
  Manager of the step by step propagation of each
  particle in a tank. In each step, determine the 
  length and generate Cherenkov photons.

*/
// ---------------------------
void ParticleInTank::DoDetSim()  
{
  if (fId == 1) GammaInt(theConfig()->SimMode);
  else if (abs(fId) == 2 || abs(fId) == 3) {
    Range=GeomLength();
    while ( DoNextStep() == CONTINUE ) {
      ComputeStepParameters();
      if ( DoProcesses() == STOP ) break;
      Propagate();
    } 
  } else return;
  
}

// ---------------------------------
/*
  Computes the distance between the
  particle and the tank.
*/
// ---------------------------------
Double_t ParticleInTank::GeomLength()
{
  Double_t aa,bb,cc,glength;
    aa = fUX*fUX+fUY*fUY;
  bb = fX*fUX+fY*fUY;
  cc = STATION_RADIUS*STATION_RADIUS-fX*fX-fY*fY;
  //   length up to the external cylinder
  if(aa!=0)glength=(-bb+sqrt(bb*bb+cc*aa))/aa;
  else glength= INFINITY;
  //   length up to the top or the bottom
  //if(fId==3)cout<<" t "<<glength<<"  "<<fX<<" "<<fY<<" "<<fZ<<" "<<fUX<<" "<<fUY<<" "<<fUZ<<" "<<aa<<" "<<bb<<" "<<cc<<endl;
  if(fUZ<0) glength = min(glength,-fZ/fUZ);
  else if(fUZ!=0)
    glength = min(glength,(STATION_HEIGHT-fZ)/fUZ);
  else  glength=2*(fX*fUX+fY*fUY);
 

  return glength;
}

// ----------------------------
/* 
   Decides if the step by step tracking must
   be continued or stopped.
*/
// ----------------------------
int ParticleInTank::DoNextStep()
{
  int ret=0;
  IStep = 0;

  switch (abs(fId)) {

  case 1:
    
    break;

  case 2:
    while (fE <= ESTEP_E[IStep]) IStep++;
    ret = (IStep<NSTEP_E && Range>=0);
    break;

  case 3:
    while (fE <= ESTEP_MU[IStep]) IStep++;
    ret = (IStep<=NSTEP_MU && Range>0);
    break;

  default:
    cerr << "Erreur dans ParticleInTank::DoNextStep..." << endl;
    // exit(0);

  }

  if (ret) ret = CONTINUE;
  else ret = STOP;
  return ret;
}

// ----------------------------------------
/*

*/
// ----------------------------------------
void ParticleInTank::ComputeStepParameters()
{
  double glength=0., energy=0.;

  switch (abs(fId)) {

  case 1:

    break;

  case 2:    
    if (IStep == 0) { // case energy > 0.1 GeV      
      DxStep=0.2;
      energy=fE;
    } else { // case energy < 0.1 GeV
      DxStep = (fE - ESTEP_E[IStep]) / DEDX_E;
      energy= (ESTEP_E[IStep-1]+ESTEP_E[IStep])/2;
    }    
    glength = GeomLength();
    Range = min(glength,DxStep);    
    Beta = sqrt(energy*(energy+2*MASS_E))/(energy+MASS_E);
    break;

  case 3:
    if (IStep==NSTEP_MU) DxStep=fE/DEDX_MU;
    else DxStep = (fE - ESTEP_MU[IStep]) / DEDX_MU;
    glength = GeomLength();
    Range = min(DxStep, glength);
    if(IStep>0) energy= (ESTEP_MU[IStep-1]+ESTEP_MU[IStep])/2;
    else energy = (10.+ESTEP_MU[IStep])/2.;
    Beta = sqrt(energy*(energy+2*MASS_MU))/(energy+MASS_MU);
    break;

  default:
    cerr << "Erreur dans ParticleInTank::ComputeStepParameters..." << endl;
    // exit(0);    

  }  
}


// ------------------------------
/*
  Simulation of processes in each step.
  For Detailed Mode:
  ------------------
  * Gamma Simulation:
  Compton Scattering
  Pair Production

  * Electron Simulation:
  Cherenkov production
  MultipleScattering
  Bremsstrahglung
  Ionization lost

  * Muon Simulation:
  Cherenkov production
  Delta Rays
  Muon Decay

*/
// ------------------------------
int ParticleInTank::DoProcesses()
{
  double probIntInOneStep;

  switch (abs(fId)) {

  case 1:
    if (GLength < IntLength) return STOP;
    if (Rand() <= ProbPair) {
      PairProduction();
      return STOP;
    } else {
      ComptonScattering();
    }
    break;

  case 2:    
    if (theConfig()->SimMode == "DETAILED") {
      if(Range>0) CherPhot(Beta,Range);
       MultipleScattering(Beta,Range);
      probIntInOneStep=gBremTable[IStep][0]*(Range*100);
      if(Rand()<probIntInOneStep) fE -= Bremstrahlung(IStep,Range,Beta);
      fE -= Range*DEDX_E;
      if (Range < DxStep || fE < EMIN_E) return STOP;
    } else { // fast simulation
      ;
    }
    break;

  case 3:    
    if (theConfig()->SimMode == "DETAILED") {
      if(Range>0) CherPhot(Beta,Range);
      if(DxStep>0) fE -= GenerateDeltaRays(Beta,Range);
      while(fE < ESTEP_MU[IStep]) {IStep++;}
      if ( MuonDecay() == 1 ) return STOP;
    } else { // fast simulation
      ;
    }    
    fE = ESTEP_MU[IStep++];
    break;

  default:
    cerr << "Erreur dans ParticleInTank::DoProcesses..." << endl;
    // exit(0);    
  }  

  return CONTINUE;
}

// ----------------------------
/*
  Propagation in space and time
  for each elementar step
*/
// ----------------------------
void ParticleInTank::Propagate()
{
  fX += Range*fUX;
  fY += Range*fUY;
  fZ += Range*fUZ;
  fT += Range/(Beta*CLIGHT);
}


// ------------------------------------------------------
// Photon processes
// ------------------------------------------------------


const Double_t kEMinPh = 0.0004;

void ParticleInTank::PairProduction()
{

}

void ParticleInTank::ComptonScattering()
{

}

// ------------------------------------------------------
// Muon processes
// ------------------------------------------------------

Double_t  ParticleInTank::GenerateDeltaRays(Double_t beta, Double_t xstep)
{
  // step in m and T in Gev
  double gamma=1/sqrt(1.-beta*beta);
  double Tmax= MASS_E* beta* beta*gamma*gamma;
  double meanNbDR=  DELTARAYCONST*xstep/beta*(1./gTMinForDeltaRays-1./Tmax);
    int nbDR=(int)Rand("POISSON",meanNbDR);
    //int nbDR=(int)meanNbDR;
  //  cout <<nbDR << endl;
  //  cout <<  DELTARAYCONST/beta*(1./gTMinForDeltaRays-1./Tmax)<<" "<<xstep<<" "<<nbDR << endl;
  double DE=0;

  //cout<<"T  "<<gTMinForDeltaRays<<" "<<Tmax<<endl;
  for(int n=0; n< nbDR; n++)
    {
      double invT=1./gTMinForDeltaRays-Rand()*(1./gTMinForDeltaRays-1./Tmax);
      double T= 1/invT;
      double costheta=T/Tmax*sqrt((Tmax*Tmax+2*MASS_E*Tmax)/
                                  (T*T+2*MASS_E*T));
     
      //  cerr<<" do delta rays"<<endl;

      // get the distance
      double dist=xstep*Rand();

      // orientation of the incident muon in the tank frame
      double UR=sqrt(fUX*fUX+fUY*fUY);

      // compute delta ray coordinates
      double dr_x = fX + dist * fUX;
      double dr_y = fY + dist * fUY;
      double dr_z = fZ + dist * fUZ;
      double time = fT + dist/(CLIGHT*beta);

      // compute the orientation in the muon frame
      double sintheta=sqrt(1-costheta*costheta);
      double phi=2*PI*Rand();
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      double ux=sintheta*cosphi;
      double uy=sintheta*sinphi;
      double uz=costheta;
      int dr_id=2;
      double dr_ux;
      double dr_uy;
      double dr_uz;
      if(UR!=0) // muon non purely vertical
        {
          dr_ux=(ux*fUZ*fUX-uy*fUY)/UR+uz*fUX;
          dr_uy=(ux*fUZ*fUY+uy*fUX)/UR+uz*fUY;
          dr_uz=-UR*ux+fUZ*uz;
        }
      else
        {
          dr_ux=ux;
          dr_uy=uy;
          if(fUZ>0) dr_uz=uz;
          else   dr_uz=-uz;
        }

      double norm= sqrt(dr_ux*dr_ux+dr_uy*dr_uy+dr_uz*dr_uz);

      double dr_E=T+MASS_E;


      double radius =sqrt(dr_x*dr_x+dr_y*dr_y);
      //     cout<<"probleme ?"<<dr_z<<" "<<radius<<endl;
      if( (dr_z>0) && (dr_z<STATION_HEIGHT)*(1-EPSILON)  && radius <STATION_RADIUS*(1-EPSILON) )

        {
          ParticleInTank* p1 = new ParticleInTank(fISta,dr_id,fId,dr_x,dr_y,dr_z,dr_ux/norm,
						  dr_uy/norm,dr_uz/norm,time,dr_E,fW);
          p1->DoDetSim("DETAILED");
          fNCherPhot += p1->fNCherPhot;
          for(int iph= 0;iph<p1->fNCherPhot;iph++){
            fTCherPhot.push_back(p1->fTCherPhot[iph]);
            fIpm.push_back(p1->fIpm[iph]);
            fNReflexions.push_back(p1->fNReflexions[iph]);
          }
          delete p1;

          DE+=T;
        }
      //else
      //cout << " Delta Ray not simulated " << endl;

    }

  return DE;

}

double FuncDistribThetaCM(double x, double rand)
{
  return M_PI*rand-x-sin(x);
}

int ParticleInTank::MuonDecay()
{
  return 0;
}


// ------------------------------------------------------
// Electron processes
// ------------------------------------------------------

Double_t  ParticleInTank::Bremstrahlung(Int_t istep,Double_t dxstep,Double_t beta){
  //cout<<" coucou"<<endl;
  double pb;
  double  EFrac, GammaKE;
  double brE ;
  pb=Rand();
  int ind = 0;
  for(int i=1;i<NDIV;i++)
    {
      if(gBremTable[istep][i]>pb)
        {ind=i; break;}
    }
  EFrac=(double)ind/(double)NDIV;
  GammaKE=fE*EFrac;

  /* Determine approx RMS gamma emission angle */
  /* (ref Rossi page 53 with q factor=1) */
  /* In fact angles are small compared to scattering so no problem */

  double mom2 = fE*(fE+2*MASS_E);
  double sig2 = pow(.0136/beta,2) / mom2 * dxstep/.361;
  double urandom=Rand();
  while(urandom==0)
    urandom=Rand();
  double theta = sqrt(-2*sig2*log(urandom));
  double costh = cos(theta);
  double sinth = sqrt(1-costh*costh);
  double phi = 2*PI * Rand();
  double cosph = cos(phi);
  double sinph = sin(phi);

  if(GammaKE>EMINBREM*dxstep)
    {
      /* Calculate the vector components of daughter direction */
      /* and rotate ref frame such that z axis lies along old ptcle track */
      //cout<<" BBREMMMMMM "<<GammaKE<<" "<<EMINBREM<<" "<<dxstep<<endl;
      double ur=sqrt(fUX*fUX+fUY*fUY);

      double vx,vy,wx,wy,wz;
      if(ur!=0)
        {
          vx = -fUY/ur;
          vy =  fUX/ur;
          wx =         -fUZ*vy;
          wy = fUZ*vx        ;
          wz = fUX*vy-fUY*vx;
        }
      else
        {
          vx=1;
          vy=0;
          wx=0;
          wy=1;
          wz=0;
        }

      // new direction

      double brUX = costh * fUX + sinth * (cosph*vx+sinph*wx);
      double brUY = costh * fUY + sinth * (cosph*vy+sinph*wy);
      double brUZ = costh * fUZ + sinth * (         sinph*wz);


      double brX=fX;
      double brY=fY;
      double brZ=fZ;
      int brId= 1;
      double brT = fT ;
      double brW = fW;
      brE = GammaKE;

      ParticleInTank* p1 = new ParticleInTank(fISta,brId,fId,brX,brY,brZ,
                                              brUX,brUY,brUZ,brT,brE,brW);
      p1->DoDetSim("DETAILED");
      fNCherPhot += p1->fNCherPhot;
      for(int iph= 0;iph<p1->fNCherPhot;iph++){
        fTCherPhot.push_back(p1->fTCherPhot[iph]);
        fIpm.push_back(p1->fIpm[iph]);
        fNReflexions.push_back(p1->fNReflexions[iph]);
      }
      delete p1;

    }
  else
    brE=0;

  return(brE);
}

void ParticleInTank::MultipleScattering(Double_t beta,Double_t xstep)
{

  double mom2 = fE*(fE+2*MASS_E);
  double sig2 = pow(.0136/beta,2) / mom2 * xstep/.361;
  double urandom=Rand();
  while(urandom==0) urandom=Rand();
  double theta = sqrt(-2*sig2*log(urandom));
  double costh = cos(theta);
  double sinth = sqrt(1-costh*costh);
  double phi = 2*PI * Rand();
  double cosph = cos(phi);
  double sinph = sin(phi);
  // construction of two unit vectors perp. to the initial direction
  // (vx,vy,0) and (wx,wy,wz)
  double ur=sqrt(fUX*fUX+fUY*fUY);

  double vx,vy,wx,wy,wz;
  if(ur!=0)
    {
      vx = -fUY/ur;
      vy =  fUX/ur;
      wx =         -fUZ*vy;
      wy = fUZ*vx        ;
      wz = fUX*vy-fUY*vx;
    }
  else
    {
      vx=1;
      vy=0;
      wx=0;
      wy=1;
      wz=0;
    }

  // new direction
  fUX = costh * fUX + sinth * (cosph*vx+sinph*wx);
  fUY = costh * fUY + sinth * (cosph*vy+sinph*wy);
  fUZ = costh * fUZ + sinth * (         sinph*wz);

}










void ParticleInTank::DoDetSim(string simmode)
{
  double glength,range,beta,dxstep,energy=1;
  double probIntInOneStep;
  range = GeomLength();

  // if(fId==3 || fId==-3) cout<<fE<<" " <<range<<endl;
  //cout<<"primaire "<<fId<<" "<<fE<<" "<<range<<endl;
  //cout<<"id = "<<fId<<endl;
  // cout<<"new particle "<<fId<<endl;
  // gamma case
  if(fId ==1) GammaInt(simmode);

  //charged particle case
  else{
    //electron case
    //    if(abs(fId)==2 && fE >= EMIN_E){
        if(abs(fId)==2){

      // determine where to begin the steps on decreasing energy
      // Low energy case: [0.1 Gev to 0.00025 GeV] gives 41 steps
      // For higher energies, istep will be 0 and a regular step
      // equivalent to 20 cm will be done
      int istep = 0;
      while(fE <= ESTEP_E[istep]) istep++;

      // in each step: determine the length and generate Cherenkov muons

      while (istep<NSTEP_E && range>=0){
        // case energy > 0.1 GeV
        if(istep==0) {
          dxstep=0.2;
          energy=fE;
        }
        // case energy < 0.1 GeV
        else{
          dxstep = (fE - ESTEP_E[istep]) / DEDX_E;
          if(istep>0)energy= (ESTEP_E[istep-1]+ESTEP_E[istep])/2;
        }

        glength = GeomLength();
        range = min(glength,dxstep);

	//        cout<<fId<<" "<<istep<<" "<<fE<<" "<<energy<<" "<<range<<" "<<glength<<endl;
        beta = sqrt(energy*(energy+2*MASS_E))/(energy+MASS_E);

        if(simmode=="DETAILED"){
          //  follows Cherenkov photons
          if(range>0)CherPhot(beta,range);


        //   //   // simulate multiple Coulomb scattering
	    MultipleScattering(beta,range);

          // simulate Bremsstraghlung
          // if the energy > 0.1 GeV the BremsProb is the same than 0.1 GeV

	  if(istep==0) probIntInOneStep=gBremTable[istep][0]*(range*100);
          else
            // low energy values
            // probIntInOneStep is in g/(cm**2)
	    probIntInOneStep=gBremTable[istep][0]*(range*100);
	  
          if(Rand()<probIntInOneStep){
            double DE=Bremstrahlung(istep,range,beta);
            //cout<<"DE= "<<DE<<endl;
            fE-=DE;
	            }
	  
	  
        }

        else   //fast simulation
          {
            Int_t n;
            Double_t avrg_pe=(38.4+44.8*TOP_FACT)*range;
            if(fW>1)  n=(Int_t)rint(fW-.5);
            else n = (Int_t)fW;

            for (Int_t i =0;i<n;i++)
              {
                Int_t npe=(Int_t)Rand("POISSON",avrg_pe);
                Double_t wgt=1.;
                if(i==1) wgt=fW-n+1;
                this->FastSim(npe,wgt);

              }
          }
        // stop if exit from water (backscattering is neglected)
        //cout<<"stop "<<range<<" "<<dxstep<<" despues de DE= "<<fE<<endl;
        if(range<dxstep || fE<0.00025) break;

        // propagate the electron
        fX += range*fUX;
        fY += range*fUY;
        fZ += range*fUZ;
        fT += range/(beta*CLIGHT);

        // energy lost
        fE -= range*DEDX_E;

        // determine where to begin the steps to start according
        // to the new energy
        istep = 0;
        while(fE <= ESTEP_E[istep]) istep++;

        //if(fE < ESTEP_E[istep] || istep!=0) istep++;
      }

    }

    // muon case
    else{
      // if (abs(fId)==3 && fE >= EMIN_MU) {
      if (abs(fId)==3) {
        //cout<<"muon "<<fE<<endl;
        // determine where to begin the steps on decreasing energy
        int istep = 0;
        while(fE <= ESTEP_MU[istep]) istep++;

        // in each step: determine the length and generate Cherenkov muons
        while (istep<=NSTEP_MU && range>0){

          if(istep==NSTEP_MU) dxstep=fE/DEDX_MU;
          else
            dxstep = (fE - ESTEP_MU[istep]) / DEDX_MU;
	 
          glength = GeomLength();
          range = min(dxstep, glength);
	  // cout<<fId<<" "<<istep<<" "<<fE<<" "<<range<<" "<<dxstep<<" "<<glength<<endl;
	  // cout<<glength<<endl;
          // cout<<fX<<" "<<fY<<" "<<fZ<<" "<<fUX<<" "<<fUY<<" "<<fUZ<<endl;
          if(istep>0)energy= (ESTEP_MU[istep-1]+ESTEP_MU[istep])/2;
          else energy = (10.+ESTEP_MU[istep])/2.;
          beta = sqrt(energy*(energy+2*MASS_MU))/(energy+MASS_MU);
          if(simmode=="DETAILED")
            {
              //  follows Cherenkov photons
              if(range>0)CherPhot(beta,range);

               //  //delta rays
                if(dxstep>0){
		  //  cout<<energy<<" ";
                double DE=GenerateDeltaRays(beta,range);
                //cout<<"delta rays "<<DE<<endl;
                fE-=DE;
              }

              while(fE < ESTEP_MU[istep]) {istep++;}

              // muon decay if stopping within geometrical range

	      double betagamma2 = (fE/MASS_MU) * (fE/MASS_MU +2);
              double stop_dist = 0.1 * MASS_MU * betagamma2 * betagamma2;
              //    cout<<fE<<" "<<range<<" "<<range-dxstep<<" "<<stop_dist<<endl;
              short flagdecay=0;

              if (range>stop_dist && range >0)
                {
                  flagdecay=1;
                  //  cout<<"muon decay"<<endl;
                  double mdX=fX+stop_dist*fUX;
                  double mdY=fY+stop_dist*fUY;
                  double mdZ=fZ+stop_dist*fUZ;
                  int mdId;
                  if(fId>0)  mdId =2;
                  else mdId=-2;
                  double randomnumber=Rand();
                  while(randomnumber==0)
                    randomnumber=Rand();
                  double delaydecay= - LIFE_MUON * log(randomnumber);
                  double mdT = fT + delaydecay;
                  double mdW = fW;
                  double mdE = MASS_MU/2 * Rand();
                  double costheta = 2*Rand()-1;
                  double sintheta = sqrt(1-costheta*costheta);
                  double phi = 2*PI * Rand();
                  double mdUX = sintheta*cos(phi);
                  double mdUY = sintheta*sin(phi);
                  double mdUZ = costheta;
                  //cout<<mdE<<" "<<mdId<<" " ;
                  //cout<<delaydecay<<endl;
                  ParticleInTank* p1 = new ParticleInTank(fISta,mdId,fId,mdX,
                                                          mdY,mdZ,mdUX,mdUY,
                                                          mdUZ,mdT,mdE,mdW);
                  p1->DoDetSim("DETAILED");
                  //cout<<p1->fNCherPhot<<endl;
                  fNCherPhot += p1->fNCherPhot;
                  for(int iph= 0;iph<p1->fNCherPhot;iph++){
                    fTCherPhot.push_back(p1->fTCherPhot[iph]);
                    fIpm.push_back(p1->fIpm[iph]);
                    fNReflexions.push_back(p1->fNReflexions[iph]);
                  }
                  delete p1;
                }

              //end of muon decay
              if(flagdecay) return;

            }

          else   //fast simulation
            {
              Int_t n;
              Double_t avrg_pe=(38.4+44.8*TOP_FACT)*range;
              if(fW>1)  n=(Int_t)rint(fW-.5);
              else n = (Int_t)fW;

              for (Int_t i =0;i<n;i++)
                {
                  Int_t npe=(Int_t)Rand("POISSON",avrg_pe);
		  Double_t wgt=1.;
                  if(i==1) wgt=fW-n+1;
                  this->FastSim(npe,wgt);

                }
            }
          fX += range*fUX;
          fY += range*fUY;
          fZ += range*fUZ;
          fT += range/(beta*CLIGHT);
          fE = ESTEP_MU[istep];
          istep++;
        }
      }
      else
        return;
    }
  }

}


void ParticleInTank::GammaInt(string simmode)
{

  const Double_t emass  = .0005;
  const Double_t emin_ph = .0004;

  Double_t xx,yy,zz,ux,uy,uz,e_ph,timr,weight;
  Double_t glength,hlength;
  Double_t e1,e2;
  Double_t ph_al,prb_pair;
  Double_t rangph,range1,range2;
  Double_t emin,smax,efin,sfin;
  Double_t c_dth,s_dth,s_dph,c_dph,dph;
  Double_t uxy,vx,vy,vz,wx,wy,wz;
  Double_t ux_f=0.,uy_f=0.,uz_f=0.,ux_el,uy_el,uz_el;
  Double_t p1,aa,bb,cc;

  Int_t active,reject;
  Int_t itab;

  //cout<<"gamma  "<<fW<<" "<<fE<<endl;;

  range1=0;
  range2=0;

  if(fE <= emin_ph) return;

  //  initial position and direction of the photon ------------------------

  e_ph=fE;
  xx=fX;
  yy=fY;
  zz=fZ;
  ux=fUX;
  uy=fUY;
  uz=fUZ;
  timr= fT;
  weight =fW;

    int iter=0;
  active=1;

  // loop over photon interactions ---------------------------------------
  // cout<<"entree dans gamma_int"<<endl;
  while(active)
  {
    //   computation of geometric length -------------
    aa = ux*ux+uy*uy;
    bb = xx*ux+yy*uy;
    cc = STATION_RADIUS * STATION_RADIUS-xx*xx-yy*yy;

    //   length up to the external cylinder
    if(aa!=0)glength=(-bb+sqrt(bb*bb+cc*aa))/aa;
    else glength=INFINITY;

    //   length up to the top or the bottom
    if(uz < 0.) hlength=-zz/uz;
    else  hlength=(STATION_HEIGHT-zz)/uz;

    //  take the smallest value
    glength=min(glength,hlength);
    if(glength<0)cerr<<"ATTENTION2"<<endl;

    //  interaction length and pair production probability : from tabulation

    /*   itab=(Int_t)(log10(e_ph)*5 + 20.5);
    if(itab <= 0) itab=0;
    if(itab > NTAB_PHOT) itab=25;
    ph_al =   TABINTLENGTH[itab-1];
    prb_pair = TABPROBPAIR[itab-1];*/

    // Changed by Denis the king of the bug !!!!!!!!!!!!!!!!!
    // To be checked !!!!!!!!!!!!!!!!! Yeeeeeeeeeeeeeeeeehhhhhhhhhhhhhhhhhhhh

    itab=(Int_t)((log10(e_ph)+ 4.1)/0.2);
    if(itab <= 0)  {ph_al=TABINTLENGTH[0];prb_pair=TABPROBPAIR[0];}
    else{
      if(itab >= (NTAB_PHOT-1)) {ph_al=TABINTLENGTH[NTAB_PHOT-1];prb_pair=TABPROBPAIR[NTAB_PHOT-1];}
      else{
        double log10Einf=-4.1+0.2*itab;
        double log10Esup=-4.1+0.2*(itab+1);
        ph_al=TABINTLENGTH[itab]+(log10(e_ph)-log10Einf)*(TABINTLENGTH[itab+1]-TABINTLENGTH[itab])/(log10Esup-log10Einf);
        prb_pair=TABPROBPAIR[itab]+(log10(e_ph)-log10Einf)*(TABPROBPAIR[itab+1]-TABPROBPAIR[itab])/(log10Esup-log10Einf);
      }
    }

    rangph = -log(Rand()) * ph_al;

    // cout<<"activity "<<glength<<" "<<rangph<<endl;
    active = glength >= rangph;

    if(active)
      {
        xx=xx+rangph*ux;
        yy=yy+rangph*uy;
        zz=zz+rangph*uz;
        timr=timr+rangph/CLIGHT;
        cc = STATION_RADIUS * STATION_RADIUS - xx*xx - yy*yy;

        iter++;

        // pair production ----------------------------------
        // final direction of electrons is the same as incident photon
        if(Rand() <= prb_pair)
          {
            //cout<<" pair production "<<endl;
            e1=(e_ph-2*emass)*Rand();
            e2=(e_ph-2*emass)-e1;

            ParticleInTank * p1=
              new ParticleInTank(fISta,2,fId,xx,yy,zz,fUX,fUY,fUZ,timr,e1,fW);
            //cout<<"e1 "<<e1<<endl;
            p1->DoDetSim(simmode);
            fNCherPhot += p1->fNCherPhot;
            for(int iph= 0;iph<p1->fNCherPhot;iph++){
              fTCherPhot.push_back(p1->fTCherPhot[iph]);
              fIpm.push_back(p1->fIpm[iph]);
              fNReflexions.push_back(p1->fNReflexions[iph]);

            }
            delete p1;

            ParticleInTank * p2=
              new ParticleInTank(fISta,-2,fId,xx,yy,zz,fUX,fUY,fUZ,timr,e2,fW);
            //cout<<"e2 "<<e2<<endl;
            p2->DoDetSim(simmode);
            fNCherPhot += p2->fNCherPhot;
            for(int iph= 0;iph<p2->fNCherPhot;iph++){
              fTCherPhot.push_back(p2->fTCherPhot[iph]);
              fIpm.push_back(p2->fIpm[iph]);
              fNReflexions.push_back(p2->fNReflexions[iph]);

            }
            delete p2;

            active=0;
          }
        // Compton scattering ------------------------------
        else
          //  Klein-Nishina energy distribution for final photon
          // (rejection method)
          {
            //cout<<" Compton scattering "<<endl;
            emin = 1/(1/e_ph + 2/emass);
            smax = emin/e_ph + e_ph/emin + pow((1+emass/e_ph-emass/emin),2) - 1;
            //   if(iter==1 && fISta==102)cout<<"compton "<<fW<<" "<<e_ph<<" "<<prb_pair<<endl;
            reject = 1;
            while(reject)
            {
              efin = emin + Rand()*(e_ph-emin);
              sfin = efin/e_ph + e_ph/efin;
              reject = Rand() > sfin/smax;
            }
             e1=e_ph - efin;

            // if needed (active photon or electron after scattering) compute final// state
            if(e1>EMIN_E || efin> emin_ph) {
              //diffusion angles of the photon
              c_dth = 1 + emass*(1/e_ph-1/efin);
              s_dth = sqrt(1 - c_dth*c_dth);
              dph=2*PI*Rand();
              c_dph=cos(dph);
              s_dph=sin(dph);
              //new direction cosines of the photon
              uxy=sqrt(ux*ux+uy*uy);
              vx=-uy/uxy;
              vy=ux/uxy;
              vz=0.;
              wx=uy*vz-uz*vy;
              wy=uz*vx-ux*vz;
              wz=ux*vy-uy*vx;
              ux_f=c_dth*ux+s_dth*(c_dph*vx+s_dph*wx);
              uy_f=c_dth*uy+s_dth*(c_dph*vy+s_dph*wy);
              uz_f=c_dth*uz+s_dth*(c_dph*vz+s_dph*wz);
            }
            //  scattered electron
            if(e1>EMIN_E) {
              p1 = sqrt(2*emass*e1+e1*e1);
              ux_el = (e_ph*ux-efin*ux_f)/p1;
              uy_el = (e_ph*uy-efin*uy_f)/p1;
              uz_el = (e_ph*uz-efin*uz_f)/p1;
              bb = xx*ux_el+yy*uy_el;
              cc = STATION_RADIUS * STATION_RADIUS - xx*xx - yy*yy;
              aa = ux_el*ux_el+uy_el*uy_el;
              glength = (-bb+sqrt(bb*bb+cc*aa))/aa;
              if(uz_el<0.)  hlength=-zz/uz_el;
              else   hlength=(STATION_HEIGHT-zz)/uz_el;

              glength = min(glength,hlength);

              ParticleInTank * p1=
                new ParticleInTank(fISta,2,fId,xx,yy,zz,ux_el,uy_el,uz_el,timr,e1,fW);
              p1->DoDetSim(simmode);
              fNCherPhot += p1->fNCherPhot;
              for(int iph= 0;iph<p1->fNCherPhot;iph++){
                fTCherPhot.push_back(p1->fTCherPhot[iph]);
                fIpm.push_back(p1->fIpm[iph]);
                fNReflexions.push_back(p1->fNReflexions[iph]);

              }
              delete p1;

            }

            //scattered photon
            active=efin >emin_ph;

            if(active)
              {
                e_ph=efin;
                ux=ux_f;
                uy=uy_f;
                uz=uz_f;
              }
          }// end of compton scattering
      }//end of active
  } //end of while

}

// ------------------------------------------------------
// Cherenkov production
// ------------------------------------------------------

void ParticleInTank::CherPhot(Double_t beta,Double_t range)
{
  // if(fId==3 || fId==-3)cerr<<"do cherenkov "<< range<<endl;
  Double_t vx,vy,vz,wx,wy,wz,uxy;
  Double_t costh_cher,sinth_cher;
  Double_t alph,calph,salph;
  Double_t cx,cy,cz;
  Double_t cx_ph,cy_ph,cz_ph,x_ph,y_ph,z_ph,r_ph;
  Double_t prb,prb_refdif,prb_spec,abslen;
  Int_t n,nphot,ipm;
  nphot=0;
  cx=fUX;
  cy=fUY;
  cz=fUZ;

  uxy=sqrt(fUX*fUX+fUY*fUY);

  if(uxy > 0.)
    {
      vx=-fUY/uxy;
      vy=fUX/uxy;
    }
  else
    {
      vx=1.;
      vy=0.;
    }

  vz=0.;
  wx=fUY*vz-fUZ*vy;
  wy=fUZ*vx-fUX*vz;
  wz=fUX*vy-fUY*vx;

  // Cherenkov effect:
  // cout<<range<<" "<<fId<<" "<<fE<<" "<<weight<<endl;
  costh_cher=1./(WATINDEX*beta);
  //  if(fId==3 || fId==-3)cout<<costh_cher<<endl;
  if(costh_cher>1)return;

  sinth_cher=sqrt(1-costh_cher*costh_cher);
  Double_t nn=(gCherFact*sinth_cher*sinth_cher*range*fW);
  n=(Int_t)nn;
  //  if(fId==3 || fId==-3)cerr<<"bilan"<<fE<<" "<<range<<" "<<nn<<endl;
  if(Rand()<nn-n) n++;

  //New Poisson distribution in the number of Cherenkov photons
  //july 2nd,2004 (as in Geant4)
 //  Int_t mu =(Int_t)nn;
//   n = (int)(Rand("POISSON",mu));

  // cout<<range<<" "<<n<<endl;
  //  cout<<"88888 "<<n<<endl;

  for(Int_t i=0; i<n;i++)
    {
      // direction of emission
      alph = 2*PI*Rand();
      calph = cos(alph);
      salph = sin(alph);
      cx_ph = costh_cher*fUX + sinth_cher * (calph*vx+salph*wx);
      cy_ph = costh_cher*fUY + sinth_cher * (calph*vy+salph*wy);
      cz_ph = costh_cher*fUZ + sinth_cher * (salph*wz);

      // point of emission
      r_ph = range*Rand();
      x_ph = fX + r_ph*cx;
      y_ph = fY + r_ph*cy;
      z_ph = fZ + r_ph*cz;

      double wav,qe;
      int iwa=0;

      //  wavelength and absorption properties at this wavelength
      prb=Rand();
      for(Int_t j=0;j<NWAVEL;j++)
        {
          //cout<<NWAVEL<<" "<<j<<" " <<  gCherPrb[j]<<endl;
          if(prb < gCherPrb[j]) {
            iwa=j;
            break;
          }
        }
      int deb=0;
      if(iwa>0) {
          wav=(WAVELMIN+(iwa+.5)*DWAVEL)*1.e9;
          prb_refdif=LINABSMAX * RELLINABS[iwa];
          abslen = WATABSLENMAX *  RELWATABSLEN[iwa];
          prb_spec = SPECFRACMAX * RELSPECFRAC[iwa];

          qe=PMTQE[iwa];
      }
      else
        {
          deb=1;
          wav=(WAVELMIN+(iwa+.5)*DWAVEL)*1.e9;
          prb_refdif=LINABSMAX * RELLINABS[iwa];
          abslen = WATABSLENMAX *  RELWATABSLEN[iwa];
          prb_spec = SPECFRACMAX * RELSPECFRAC[iwa];
          //wav = (WAVELMIN+(j+.5)*DWAVEL)*1.e9;
          qe=PMTQE[iwa];

        }

      //cout<<iwa<<" "<<wav<<" "<<abslen<<" " <<prb_refdif<<" "<<qe<<" "<<deb<<endl;
      //   cout<<"pos "<<x_ph<<" " <<y_ph<<" "<<z_ph<<" " <<cx_ph<<" " <<cy_ph<<" "<<cz_ph<<endl;

      // following a Cherenkov photon
      ipm=this->Follow(x_ph,y_ph,z_ph,cx_ph,cy_ph,cz_ph,
                       prb_refdif,prb_spec,abslen);
      if(ipm>0) { nphot++;}

    }

  fNCherPhot+= nphot;

  //cout<<" "<<nphot<<endl;
  // if(fId==3 || fId==-3)cout <<fId<<" "<< n << " "<<nphot <<" "<<weight<<" "<<e1<<" "<<range<<endl;
  //if(fId==3 || fId==-3) cout <<fISta<<" "<<fId<<" "<< n << " "<<nphot <<" "<<fNCherPhot<<" "<<range<<" "<<beta<<" "<<sinth_cher<<endl;
  // if(fId==3 || fId==-3) cout<<"cer "<<n<<" "<<nphot<<" "<<fNCherPhot<<" "<<endl;
  // if(fId==3)cout<<endl;
}



Int_t ParticleInTank::Follow(Double_t x,Double_t y,Double_t z,
                           Double_t cosx,Double_t cosy,Double_t cosz,
                           Double_t prb_refdif,Double_t prb_spec,
                           Double_t abslen)
{
  Double_t zc_pm;
  Double_t xx,yy,zz,cxx,cyy,czz,norm,prod;
  Double_t aa,bb,cc,glength,zs;
  Double_t c2bet,cbet,sbet,alph,calph,salph;
  Int_t active;
  Double_t zhit;
  Int_t ipm;
  Double_t tlength;

  zc_pm=STATION_HEIGHT+HC_PM;

  ipm=0;
  tlength=0.;
  xx=x;
  yy=y;
  zz=z;
  cxx=cosx;
  cyy=cosy;
  czz=cosz;
  norm=sqrt(cxx*cxx+cyy*cyy+czz*czz);
  cxx=cxx/norm;
  cyy=cyy/norm;
  czz=czz/norm;
  // prb_spec=0.1;
  active=1;
  Int_t ind;
  ind=0;

 // computation of geometric length
  while(active)
    {

      aa=cxx*cxx+cyy*cyy;
      bb=xx*cxx+yy*cyy;
      cc=STATION_RADIUS*STATION_RADIUS-xx*xx-yy*yy;
      glength=(-bb+sqrt(bb*bb+aa*cc))/aa;
      zs=zz+glength*czz;

      //   side wall  -------------------------------------
      if(zs > 0.&& zs < STATION_HEIGHT) {

        // absorption in water or in wall
        active = Rand() < exp(-glength/abslen);
        active = active && Rand()< prb_refdif;

        if(active) {
          xx += glength*cxx;
          yy += glength*cyy;
          zz = zs;
          tlength += glength;
          ind++;
          //  Double_t intheta= acos(STATION_HEIGHT/tlength)*RAD2DEG;
          /*   geometry of reflection/diffusion unit vectors (w.r.t. surface)
               normal  n  : -xx/STATION_RADIUS  -yy/STATION_RADIUS   0
               tangent t1 :  yy/STATION_RADIUS  -xx/STATION_RADIUS   0
               t2 :         0               0                1  */

          // specular reflection
          if(Rand() < prb_spec)
            {
              prod = (xx*cxx+yy*cyy)/STATION_RADIUS;
              cxx += -2*prod*xx/STATION_RADIUS;
              cyy += -2*prod*yy/STATION_RADIUS;


              // diffusion (sinus.cosinus law -->  uniform distr. in cos(2*th)
            }
          else
            {
              c2bet = 2*Rand()-1;
              cbet = sqrt((1+c2bet)/2);
              sbet = sqrt((1-c2bet)/2);
              alph = 2*Rand()*PI;
              calph = cos(alph);
              salph = sin(alph);
              cxx = (-cbet*xx+sbet*salph*yy)/STATION_RADIUS;
              cyy = (-cbet*yy-sbet*salph*xx)/STATION_RADIUS;
              czz = sbet * calph;
              //cout<<"1"<<" "<<intheta<<" "<<180-acos(cxx)*RAD2DEG<<endl;

            }
          // if(ind>=1)cout<<ind<<" side "<< czz<<endl;
          // cout<<ind<<" side "<< czz<<endl;
        }//end active
      }//end side wall

      //  top wall ---------------------------------------------------
      else
        if(czz > 0.) {
          glength = (STATION_HEIGHT-zz)/czz;
          active = Rand() < exp(-glength/abslen);

          if(active) {
            xx += glength*cxx;
            yy += glength*cyy;
            zz = STATION_HEIGHT;
            tlength += glength;


            //  PMT simulation portion of a semisphere PMT
            //  PMT entry ?
            for(Int_t j=0;j<NPM;j++){
              bb = cxx*(xx-X_PM[j])+cyy*(yy-Y_PM[j])+czz*(zz-zc_pm);
              cc = RAD_PM * RAD_PM - pow((xx-X_PM[j]),2) - pow((yy-Y_PM[j]),2)
                -pow((zz-zc_pm),2);

              if(bb*bb+cc >0.)
                {
                  zhit = zz-(bb+sqrt(bb*bb +cc))*czz;
                  if(zhit <= STATION_HEIGHT &&
                     zhit >( STATION_HEIGHT + HC_PM-RAD_PM))
                    {
                      ipm=j+1;
                      fTCherPhot.push_back(tlength*WATINDEX/CLIGHT+fT);
                      fIpm.push_back(ipm);
                      fNReflexions.push_back(ind);
                      // cout<<ipm<<" "<<ind<<" "<<tlength<<endl;
                      // cout<<ind<<endl;
                      return(ipm);

                    }
                }
            }//loop on PMT

            // cout<<"999999 "<<ipm<<" "<<ind<<" "<<tlength<<" "<<xx<<" "<<yy<<endl;
            ind++;

            // no PMT entry
            active=Rand() < TOP_FACT * prb_refdif;
            if(active){
              // specular reflection

              if(Rand() < prb_spec)
                czz=-czz;
              else
                // diffusion (sinus.cosinus law -->  uniform distr. in cos(2*th)
                {
                  c2bet = 2*Rand()-1;
                  sbet = sqrt((1-c2bet)/2);
                  alph = Rand() * 2 * PI;
                  cxx = sbet * cos(alph);
                  cyy = sbet * sin(alph);
                  czz = -sqrt((1+c2bet)/2);
                }
              //if(ind>=1)cout<<ind<<" top "<< czz<<" " <<ipm<<endl;

            }//end active

          }//end active
	  
        }//end top
      
      // bottom wall ----------------------------------------------------
      
	else {
	  
          glength = -zz/czz;
          active = Rand()< exp(-glength/abslen);
          active = active && Rand() < prb_refdif;
          // Double_t thetain = acos(-czz)*RAD2DEG;
          if(active){
            xx += glength*cxx;
            yy += glength*cyy;
            zz=0.;
            tlength += glength;
            ind++;
            Double_t spec = Rand();
            // Double_t thetain = acos(-czz)*RAD2DEG;
            if(spec < prb_spec){
              czz=-czz;
              //  cout<<"0 "<<thetain<<" "<<acos(czz)*RAD2DEG<<endl;
            }
            else
              // diffusion (sinus.cosinus law -->  uniform distr. in cos(2*th))
              {
                c2bet = 2*Rand()-1;
                sbet = sqrt((1-c2bet)/2);
                alph = Rand()*2*PI;
                cxx = sbet * cos(alph);
                cyy = sbet * sin(alph);
                czz = sqrt((1+c2bet)/2);

                //cout<<"1 "<<thetain<<" "<< acos(czz)*RAD2DEG<<endl;
                //cout<<c2bet<<" "<<sbet<<" "<<alph<<" "<<cxx<<" "<<cyy<<" "<<czz<<endl;
              }

            //  if(ind>=1)cout<<ind<<" bottom "<< czz<<endl;
	    // cout<<"999999 "<<ipm<<" "<<ind<<" "<<tlength<<" "<<xx<<" "<<yy<<endl;

          }//end active

        }//end bottom */
    }
  // cout<<ind<<" "<<tlength<<" "<<ipm<<endl;
  return(ipm);
}



void ParticleInTank::FastSim(Int_t npe,Double_t wgt)
{

  Int_t npeeff,djt,irb;
  Double_t weff,prb;
  Int_t i,j;
  Double_t bincontent[MAXNUMBEROFADCBINS];


  for(Int_t i = 0;i<MAXNUMBEROFADCBINS;i++)
    bincontent[i]=0;


  if (npe <= 0) return;

  npeeff = (Int_t) ((Double_t)npe * gNpesamples);
  weff   = wgt / gNpesamples;
  djt    = (Int_t) (fT /SAMPLINGTIME);


  for ( i = 0;i<npeeff;i++)
    {
      prb = Rand();
      irb = (Int_t) (gNshapezones * prb);


      for ( j = gShapezone[irb] ; j<gNshape;j++)
          if(prb<gIntshape[j]) break;


      j+= djt;



      if (j > 0 && j<= gNbins)
        {

          bincontent[j]=weff+ bincontent[j];

        }

    }


  if(fBinContent.size()==(unsigned)gNbins){
   for ( j = 0 ; j<gNbins ; j++)
     {
       fBinContent[j]=(bincontent[j]+weff*(0.5-Rand()));

     }
  }
  else if (fBinContent.size()==0)
    {
      for ( j = 0 ; j<gNbins ; j++)
        {
          fBinContent.push_back(bincontent[j]+weff*(0.5-Rand()));
        }
    }
  else
    cerr<<" Pb of vector size fBinContent in  ParticleInTank"<<endl;

}
