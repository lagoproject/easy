//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File BuildProcessesTables.h 

 */
//----------------------------------------------------------------------

#ifndef BUILDPROCESSESTABLES_H
#define BUILDPROCESSESTABLES_H

//CONSTANTS RELATED TO PHOTON INTERACTION IN TANK

const int NTAB_PHOT = 25;

const double TABINTLENGTH[NTAB_PHOT] =
{.06,.07,.08,.10,.12,.18,.20,.27,.30,.40,.50,.55,.59,.60,.59,
 .56,.55,.52,.51,.50,.48,.47,.46,.46,.46};

const double TABPROBPAIR[NTAB_PHOT] =
{ 0., 0., 0., 0., 0., 0., 0.,.03,.08,.17,.33,.42,.59,.70,.79,
  .88,.91,.93,.95,.97,.99, 1., 1., 1., 1.};

//Tables size for Bremstrahlung

const int  NDIV                  = 100;
const double EMINBREM            = 0.100;  //Energy minimum for bremstrahlung in GeV/m

// Table for electron energy distribution for muon decay
const int NENER_ELEC_MUDECAY = 100;



//ENERGY BINS FOR ELECTRONS and MUONS ( in GeV)

const int  NSTEP_MU           = 43;
const int  NSTEP_E            = 41;

const double ESTEP_MU[NSTEP_MU] = {
  10.,9.,8.,7.,6.,5.,4.,3.,2.,
  1.,.95,.9,.85,.8,.75,.7,.65,.6,
  .55,.5,.45,.4,.35,.3,.25,.2,.18,
  .16,.14,.12,.1,.09,.08,.07,.06,.05,
  0.04,0.03,0.02,0.01,0.005,0.001,0.0005};

const double ESTEP_E[NSTEP_E] = {
  0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,
  .01,.0095,.009,.0085,.008,.0075,.007,.0065,
  .006,.0055,.005,.0045,.004,.0035,.003,.0025,
  .002,.0018,.0016,.0014,.0012,.001,.0009,.0008,
  .0007,.0006,.0005,.00045,.0004,.00035,.0003,.00025};

const double WATINDEX            = 1.33;
const double WATABSLENMAX        = 30;
const double SPECFRACMAX         = 1.;
const double LINABSMAX           = 0.973;


const int NWAVEL                 = 30;
const double WAVELMIN            = 295.E-9;
const double DWAVEL              = 10.E-9;

//Specular reflexion limit


const double RELWATABSLEN[NWAVEL]  = {
  0.138054 , 0.171911 , 0.21003 , 0.256239 , 0.308537 ,
  0.37735 , 0.484984 , 0.659117 , 0.753154 , 0.842399 , 
0.888049 , 0.94298 , 1.00022 , 0.993223 , 0.989593 , 
0.928307 , 0.918782 , 0.819461 , 0.727707 , 0.578173 , 
0.390255 , 0.300568 , 0.284079 , 0.256601 , 0.229324 ,
 0.205447 , 0.183962 , 0.12614 , 0.0299486 , 0 }; 


const float RELSPECFRAC[NWAVEL]  =

{0.20227 , 0.203915 , 0.204977 , 0.20555 , 0.205714 , 0.205533 , 0.205063 , 0.20435 , 0.203434 , 0.202347 , 0.201118 , 0.199771 , 0.198327 , 0.196803 , 0.195214 , 0.193573 , 0.191891 , 0.190178 , 0.188441 , 0.186689 , 0.184925 , 0.183157 , 0.181388 , 0.179622 , 0.177862 , 0.176111 , 0.174372 , 0.172646 , 0.170935 , 0.169241 };


const double RELLINABS[NWAVEL]     = 
{.7755,.8178,.8541,.8849,.9107,
 .9333,.9500,.9642,.9753,.9840,
 .9903,.9948,.9976,.9993,1.,
 1.,.9993,.9985,.9975,.9963,
 .9953,.9945,.994,.9937,.9937,
 .9939,.9942,.9947,.9953,.9957};

const double COLLEF              = 0.9;

const double PMTQE[NWAVEL]      = {
  .04, .09, .14, .17, .21, .22, .23, .24, .24, .24,
  .23, .225, .22, .21, .20, .19, .18, .16, .15, .13,
  .11, .10, .09, .08, .07, .06, .05, .04, .03, .02};


void BuildProcessesTables();
double CalcPhiUNudNu(double U, double Nu, double dNu, double A, double Z);

#endif
