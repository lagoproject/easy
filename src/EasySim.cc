//----------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  

 */
//----------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <pthread.h>
#include "EasySim.h"
#include "Array.h"
#include "Station.h"
#include "Event.h"
#include "Calib.h"
#include "Utils.h"
#include "ManageEvent.h"
#include "BuildProcessesTables.h"
#include "ShowerParam.h"
#include "Constants.h"


//----------------------------------------------------------------
/*
  global variables
*/
//----------------------------------------------------------------
TROOT root ("Main", "Tank Simulation");	// ROOT workspace
//Array *gCoreArray;
HitStation *gHStation;
extern double gMuonEnergyAngularSpectrum[ANGULAR_FLUX_SIZE];
extern double gMuonProbaAngularSpectrum[ANGULAR_FLUX_SIZE];

void new_read_config() {
  ifstream cfgfile;
  cfgfile.open ("./default.inp", ios::in);
  theConfig()->Mode="CALIB";
  theConfig()->SimMode="DETAILED";
  theConfig()->MuEmMode="DEFAULT";
  theConfig()->ElecMode="PM";
  theConfig()->NEvents=100;
  theConfig()->RunNumber="01";
  theConfig()->PartMode="RANDOM";
  theConfig()->PartCode=0;
  theConfig()->PartEnergy=0;
  theConfig()->PartTheta=0;
  theConfig()->PartAddMode="DEFAULT";
  theConfig()->NPartMultiple=20;
  char keyword[1024];
  while(cfgfile >> keyword) {
    if (keyword[0]=='#') cfgfile.ignore (1000, '\n');
    else if (strcmp(keyword,"SimMode")==0) cfgfile >> theConfig()->SimMode;
    else if (strcmp(keyword,"MuEmMode")==0) cfgfile >> theConfig()->MuEmMode;
    else if (strcmp(keyword,"ElecMode")==0) cfgfile >> theConfig()->ElecMode;
    else if (strcmp(keyword,"NEvents")==0) cfgfile >> theConfig()->NEvents;
    else if (strcmp(keyword,"RunNumber")==0) cfgfile >> theConfig()->RunNumber;
    else if (strcmp(keyword,"PartMode")==0) cfgfile >> theConfig()->PartMode;
    else if (strcmp(keyword,"PartCode")==0) cfgfile >> theConfig()->PartCode;
    else if (strcmp(keyword,"PartEnergy")==0) cfgfile >> theConfig()->PartEnergy;
    else if (strcmp(keyword,"PartTheta")==0) cfgfile >> theConfig()->PartTheta;
    else if (strcmp(keyword,"PartAddMode")==0) cfgfile >> theConfig()->PartAddMode;
    else if (strcmp(keyword,"NPartMultiple")==0) cfgfile >> theConfig()->NPartMultiple;
  }
  cfgfile.close();
}

//----------------------------------------------------------------
/*
  reading config file
*/
void
read_config ()
{
  ifstream cfgfile;
  cfgfile.open ("./default.inp", ios::in);
  // Reading the parameters from default.inp file
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->Mode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->SimMode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->MuEmMode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->ElecMode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->NEvents;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->RunNumber;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->PartMode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->PartCode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->PartEnergy;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->PartTheta;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->PartAddMode;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile >> theConfig ()->NPartMultiple;
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');
  cfgfile.ignore (100, '\n');

  if (theConfig ()->Mode == "SHOWER") {
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->ArrayMode;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->ArrayFileName;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->EventMode;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->Acceptance;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->EasAvrg;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->EasSpread;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->NorAvrg;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->NorSpread;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->PhiRotation;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->RMin;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->RMax;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->DeltaR;
    cfgfile.ignore (100, '\n');
    cfgfile.ignore (100, '\n');
    cfgfile >> theConfig ()->DeltaAzim;

  }

  theConfig ()->DeltaAzim = theConfig ()->DeltaAzim * DEG2RAD;

  cfgfile.close ();
}


//----------------------------------------------------------------
/*
  parsing arguments
*/
void
parse_argv (int argc, char **argv)
{
  theConfig ()->NoFileName = theConfig ()->TwoFileName =
    theConfig ()->AsciiFile = 0;
  switch (argc) {
  case 1:
#ifndef CALIBONLY
    cerr << endl << "No file entered: OK IF CALIBRATION MODE" << endl << endl;
#endif
    theConfig ()->NoFileName = 1;
    break;
  case 2:
    cerr << "\tSelected input file : " << endl << "\t" << argv[1] << endl;
    theConfig ()->InRootFileName = argv[1];
    break;
  case 3:
    cerr << "\tSelected output file: " << endl << "\t" << argv[2] << endl;
    theConfig ()->InRootFileName = argv[1];
    theConfig ()->OutRootFileName = argv[2];
    theConfig ()->TwoFileName = 1;
    break;
  case 4:
    cerr << "\tSelected output file: " << endl << "\t" << argv[2] << endl;
    theConfig ()->InRootFileName = argv[1];
    theConfig ()->OutRootFileName = argv[2];
    theConfig ()->TwoFileName = 1;
    theConfig ()->AsciiFileName = argv[3];
    theConfig ()->AsciiFile = 1;
    break;
  default:
    cerr << " ++++ Too many arguments ++++" << endl;
    cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cerr << "Usage : " << argv[1] << " Grnd-Root-filename " << endl;
    cerr <<
      " \t where Grnt-Root-filename is the ROOT file containing the Ground particles  "
      << endl;
    cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    exit (1);
  }

}


//----------------------------------------------------------------
/*
  Building out file names, ...
*/
void
build_outputs (Event * event)
{
  if (theConfig ()->Mode == "CALIB") {
    for (int ii = 0; ii < 80; ii++)
      cerr << "-";
    cerr << endl;
    cerr << "\t\tCALIBRATION MODE " << endl;
    theConfig ()->OutRootFileName = "calib";
    if (theConfig ()->PartMode == "VEM")
      theConfig ()->OutRootFileName += "_vem";
    if (theConfig ()->PartMode == "RANDOM")
      theConfig ()->OutRootFileName += "_rand";
    if (theConfig ()->PartMode == "SCINTILLATOR")
      theConfig ()->OutRootFileName += "_scint";
    if (theConfig ()->PartAddMode == "MULTIPLE")
      theConfig ()->OutRootFileName += "_mul";
    if (theConfig ()->PartAddMode == "DOUBLE")
      theConfig ()->OutRootFileName += "_double";
    if (theConfig ()->PartMode == "HORIZONTAL")
      theConfig ()->OutRootFileName += "_hori";
    if (theConfig ()->PartMode == "FIXEDTHETA")
      theConfig ()->OutRootFileName += "_theta";
    if (theConfig ()->SimMode == "FAST")
      theConfig ()->OutRootFileName += "_fast";
    if (theConfig ()->SimMode == "DETAILED")
      theConfig ()->OutRootFileName += "_det";
    theConfig ()->OutRootFileName =
      theConfig ()->OutRootFileName + "_" + theConfig ()->RunNumber;
    theConfig ()->OutRootFileName = theConfig ()->OutRootFileName + ".root";
  } else if (theConfig ()->Mode == "SHOWER") {
    for (int ii = 0; ii < 80; ii++)
      cerr << "-";
    cerr << endl;
    cerr << "\t\tSHOWER MODE " << endl;
    if (!theConfig ()->TwoFileName) {
      if (theConfig ()->NoFileName == 1) {
	cerr <<
	  "Please enter an ouput filename or change the mode of simulation" <<
	  endl;
	exit (0);
      }
      theConfig ()->OutRootFileName = theConfig ()->InRootFileName;

      while (string_contains (theConfig ()->OutRootFileName, "/")) {
	theConfig ()->OutRootFileName =
	  string_after (theConfig ()->OutRootFileName, "/");
      }

      if (theConfig ()->SimMode == "FAST")
	theConfig ()->OutRootFileName =
	  "Sim_fast_" + theConfig ()->OutRootFileName;
      else if (theConfig ()->SimMode == "SAMPLE")
	theConfig ()->OutRootFileName =
	  "Sample_" + theConfig ()->OutRootFileName;
      else
	theConfig ()->OutRootFileName =
	  "Sim_" + theConfig ()->OutRootFileName;

    }
  } else {
    cerr << "Please choose a mode of simulation (CALIB | SHOWER)" << endl;
    exit (0);
  }

  theConfig ()->OutRootFile =
    new TFile (theConfig ()->OutRootFileName.c_str (), "RECREATE");

  if (theConfig ()->Mode == "SHOWER"
      && theConfig ()->ArrayMode == "STARARRAY") {
    theConfig ()->OutTree =
      new TTree ("TreeEvent", "tree for HitStation objects");
    gHStation = new HitStation;
    theConfig ()->OutBranch =
      theConfig ()->OutTree->Branch ("hstation", "HitStation", &gHStation,
				     16000, 0);
    delete gHStation;
    //      outltp.open("ltp.dat",ios::out);
  } else {
    theConfig ()->OutTree = new TTree ("TreeEvent", "tree for Event objects");
    event = new Event;
    theConfig ()->OutBranch =
      theConfig ()->OutTree->Branch ("event", "Event", &event, 16000, 0);
    delete event;
  }


  //  cerr<<endl;for (int ii=0; ii<80; ii++) cerr<<"-";cerr<<endl;

}

//----------------------------------------------------------------
/*
  Building out file names, ...
*/
int
CheckHexaCore (Event * event)
{
#define HEXADIST(i) sqrt(pow(gHexagArray->fStationList[i].fNorthing-event->fNorCore,2.)+pow(gHexagArray->fStationList[i].fEasting-event->fEasCore,2.))
  double distref =
    sqrt (pow (theConfig ()->NorAvrg - event->fNorCore, 2.) +
	  pow (theConfig ()->EasAvrg - event->fEasCore, 2.));
  for (unsigned int i = 0; i < gHexagArray->fStationList.size (); i++)
    if (HEXADIST (i) < distref) {
      cerr << "\tDiscarding: " << distref << "  " << HEXADIST (i) << endl;
      return 0;
    }

  return 1;
}


void toto(void* a)
{
	a=a;
}
//----------------------------------------------------------------
/*
  main function
*/
int
main (int argc, char **argv)
{
  Event *event = NULL;
  int nt3 = 0, nt4 = 0;
  // initialisations
  cerr << endl;
  for (int ii = 0; ii < 80; ii++)
    cerr << "-";
  cerr << endl << endl;
#ifndef CALIBONLY
  cerr << "\t\t\t\tPROGRAM LAGO EASY" << endl;
  cerr << "\t\t\t\t+-+-+-+-+-+-+-+-+" << endl;
  cerr << "\t\tComplete simulation of the LAGO detector" << endl;
#else
  cerr << "\t\t\tLAGO EASY" << endl;
  cerr << "\t\t\tBased on Auger EASYSIM" << endl;
  cerr << "\t\t\t+-+-+-+-+-+-+-+-" << endl;
  cerr << "\t\tSimplified simulation of the LAGO detector" << endl;
#endif
  float vol=pow( STATION_RADIUS,2) * PI * STATION_HEIGHT;
  float surf= RAD_PM *( RAD_PM-HC_PM) *2 * PI*10000.;
  cerr << "\t\tVolume of each tank = " <<vol<<"  m3"<< endl;
  cerr << "\t\tEach tank has = " <<NPM<<"  PMT"<< endl;
  cerr << "\t\tSurface of each PMT = " <<surf <<"  cm2 " <<endl;
  cerr << endl;
  for (int ii = 0; ii < 80; ii++)
    cerr << "-";
  cerr << endl << endl;
#ifndef CALIBONLY
  read_config ();
#else
  new_read_config();
#endif
  parse_argv (argc, argv);
  build_outputs (event);

  // building some tables
  if (theConfig ()->SimMode == "DETAILED") {
    cerr << "\tDETAILED SIMULATION CHOSEN" << endl << endl;
    BuildProcessesTables ();
  } else if (theConfig ()->SimMode == "FAST") {
    cerr << "\tFAST SIMULATION CHOSEN" << endl;
    //BuildIntegratedMuonPulse(theConfig()->Mode);
  } else if (theConfig ()->Mode == "SHOWER"
	     && theConfig ()->SimMode == "SAMPLE") {
    cerr << "\tSIMULATION WILL STOP AFTER SAMPLING" << endl;
  } else {
    cerr << " You have not chosen an existing mode of simulation" << endl;
    exit (0);
  }
  for (int ii = 0; ii < 80; ii++)
    cerr << "-";
  cerr << endl;

  // mode shower
  if (theConfig ()->Mode == "SHOWER") {
    // read simulated shower
    gShParamP = new ShowerParam (theConfig ()->InRootFileName);

    // define array
    if (theConfig ()->ArrayMode == "READFILE") {
      gArray = new Array (theConfig ()->ArrayFileName);
    } else if (theConfig ()->Acceptance) {
      gArray = new Array (atoi (theConfig ()->ArrayFileName.c_str ()));
      vector < Station > hexa;
      for (unsigned int i = 0; i < gArray->fStationList.size (); i++) {
	if (gArray->fStationList[i].fId == 16 ||
	    gArray->fStationList[i].fId == 17 ||
	    gArray->fStationList[i].fId == 21 ||
	    gArray->fStationList[i].fId == 23 ||
	    gArray->fStationList[i].fId == 26 ||
	    gArray->fStationList[i].fId == 27) {
	  hexa.push_back (gArray->fStationList[i]);
	}
      }
      gHexagArray = new Array (hexa);
    } else if (theConfig ()->ArrayMode == "STARARRAY") {
      gArray =
	new Array (gShParamP->Theta, gShParamP->Phi, theConfig ()->EasAvrg,
		   theConfig ()->NorAvrg);
    } else {
      gArray = new Array (atoi (theConfig ()->ArrayFileName.c_str ()));
    }


    //average position and spreading limits of the events
    if (theConfig ()->EventMode == "DEFAULT") {
      theConfig ()->EasAvrg = gArray->fEasMean;
      theConfig ()->NorAvrg = gArray->fNorMean;
      theConfig ()->EasSpread *= gArray->fEasWidth / 2.;
      theConfig ()->NorSpread *= gArray->fNorWidth / 2.;
    }

    for (Int_t n = 0; n < theConfig ()->NEvents;) {
      event = GenerateEvent (n);
      if (theConfig ()->Acceptance && CheckHexaCore (event) == 0)
	continue;
      n++;
      gShParamP->XCore = event->fXCore;
      gShParamP->YCore = event->fYCore;
      gShParamP->EasCore = event->fEasCore;
      gShParamP->NorCore = event->fNorCore;

      vector < SampStation > samp_list;
      DoSampling (event->fEasCore, event->fNorCore, &samp_list);
      if (theConfig ()->ArrayMode != "STARARRAY") AddRandomParticles (&samp_list);
      UInt_t samp_size = samp_list.size ();
      for (UInt_t ist = 0; ist < samp_size; ist++) {
	HitStation *hst = GenerateHitStation (samp_list[ist]);
	if (!hst)
	  continue;
	if (theConfig ()->ArrayMode == "STARARRAY") {
	  //write in ltp.dat if it is a tot
	  //              if((hst->fT2ToT)==2)outltp<<hst->fId<<" ";
//	  theConfig ()->OutBranch->SetAddress (&hst);
	  gHStation = hst;
	  theConfig ()->OutTree->Fill ();
	} else {
	  event->AddTank (hst);
	}

	delete hst;

      }

      if (theConfig ()->ArrayMode != "STARARRAY") {
	event->CountTanks ();
	nt3 += event->DoT3Sim ();
	nt4 += event->DoT4Sim ();
	theConfig ()->OutRootFile->cd ();
	theConfig ()->OutBranch->SetAddress (&event);
	theConfig ()->OutTree->Fill ();
	cerr << endl << "\t---> Writing event in root output file" << endl;
      }

      delete event;
      cerr << endl;
      for (int ii = 0; ii < 80; ii++)
	cerr << "-";
      cerr << endl;
    }


    delete gShParamP;
  }
  // mode calib
  else if (theConfig ()->Mode == "CALIB") {
    Int_t id = theConfig ()->PartCode;
    Double_t en, w;
    Double_t theta = 0, azim = 0, thetamoy = 0;
    Double_t time, time0;
    double d_perp, d_rand;
    Double_t xx, yy, zz, cx, cy, cz, uxy2;
    Double_t x, y;
    Double_t eps;
    double prb;
    int nside;
    int ntop;
    int ebin;
    int npart;
    ParticleInTank *part;

    HitStation *hst;

    eps = 1.e-3;
    nside = 0;
    ntop = 0;
    time0 = 0;
    time = 0;
    ebin = 0;


    //* Calibration with random muons *//
    if (theConfig ()->PartAddMode == "MULTIPLE")
      npart = theConfig ()->NPartMultiple;
    else if (theConfig ()->PartAddMode == "DOUBLE")
      npart = 2;
    else
      npart = 1;

    cerr << "Simulation of " << theConfig ()->NEvents;
    if (theConfig ()->PartCode == 0)
      cerr << "  muons, electrons, photons ";
    else if (theConfig ()->PartCode == 1)
      cerr << " photons ";
    else if (theConfig ()->PartCode == 2)
      cerr << " electrons ";
    else if (theConfig ()->PartCode == 3)
      cerr << " muons ";
    else {
      cerr << "Unknown type of particles" << endl;
      exit (0);
    }

    if (theConfig ()->PartEnergy == 0)
      cerr << " spectrum" << endl;
    else
      cerr << " monoenergetic " << theConfig ()->PartEnergy << " GeV" << endl;
    if (theConfig ()->PartMode == "VEM")
      cerr << "central and vertical" << endl;
    else if (theConfig ()->PartMode == "SCINTILLATOR")
      cerr << "according to scintillators positions" << endl;
    else if (theConfig ()->PartMode == "FIXEDTHETA")
      cerr << "at fixed zenital angle " << theConfig ()->
	PartTheta << " degrees and all around the tank" << endl;
    else if (theConfig ()->PartMode == "HORIZONTAL")
      cerr << "horizontal around the tank" << endl;
    else if (theConfig ()->PartMode == "RANDOM")
      cerr << "randomly in sin2theta" << endl;
    else {
      cerr << "Unknown mode" << endl;
      exit (0);
    }

    if (theConfig ()->PartAddMode == "DOUBLE")
      cerr << "1 additionnal particle delayed per event" << endl;
    if (theConfig ()->PartAddMode == "MULTIPLE")
      cerr << theConfig ()->
	NPartMultiple << " simultanous particles per event" << endl;

    cerr << "----------------------------------------------------" << endl;
    for (Int_t i = 0; i < theConfig ()->NEvents; i++) {
      for (Int_t ipart = 0; ipart < npart; ipart++) {
	ebin = 0;
	prb = 1;
	w = 1;
	
	
	if (theConfig ()->PartEnergy > 0)
	  en = theConfig ()->PartEnergy;
	else			// Generation of particles with different energies
	  if (theConfig ()->PartCode > 0) {	//one kind of particle
	    if ((theConfig ()->PartMode != "RANDOM")
		&& (theConfig ()->PartCode == 3)) {
	      prb = 1;
	      //while(prb>0.498)
	      prb = Rand ();
	      for (Int_t j = 0; j < ANGULAR_FLUX_SIZE; j++) {
		if (prb < gMuonProbaAngularSpectrum[j]) {
		  ebin = j;
		  break;
		}
	      }
	      if (ebin == 0)
		en = gMuonEnergyAngularSpectrum[ebin];
	      else {
		en = ((gMuonProbaAngularSpectrum[ebin - 1]) *
		      gMuonEnergyAngularSpectrum[ebin]
		      - (gMuonProbaAngularSpectrum[ebin]) *
		      gMuonEnergyAngularSpectrum[ebin - 1]
		      + (prb) * (gMuonEnergyAngularSpectrum[ebin - 1] -
				 gMuonEnergyAngularSpectrum[ebin]))
		  / ((gMuonProbaAngularSpectrum[ebin - 1]) -
		     (gMuonProbaAngularSpectrum[ebin]));
	      }
	      cout<<ebin<<" "<<en<<endl;
	    } else {
	      prb = 1;
	      //while(prb>0.498)
	      prb = Rand ();
	      for (Int_t j = 0; j < FLUX_SIZE; j++) {
		if (prb < FLUX_INT[id - 1][j]) {
		  ebin = j;
		  break;
		}
	      }
	      if (ebin == 0)
		en = FLUX_EN[id - 1][ebin];
	      else {
		en = ((FLUX_INT[id - 1][ebin - 1]) * FLUX_EN[id - 1][ebin]
		      - (FLUX_INT[id - 1][ebin]) * FLUX_EN[id - 1][ebin - 1]
		    + (prb) * (FLUX_EN[id - 1][ebin - 1] -
			       FLUX_EN[id - 1][ebin]))
		  / ((FLUX_INT[id - 1][ebin - 1]) - (FLUX_INT[id - 1][ebin]));
	      }
	    }
	    //  cout<<en<<endl;
	      }
	  else {		// mu,el,photons in the flux
	    double probpart = Rand ();
	    if (probpart < FLUX_RATIO[0])
	      id = 1;		//photons
	    else if (probpart < FLUX_RATIO[1])
	      id = 2;		//electrons
	    else
	      id = 3;		//muons
	    prb = Rand ();
	    
	    ebin = 0;
	    for (Int_t j = 0; j < FLUX_SIZE; j++) {
	      if (prb < FLUX_INT[id - 1][j]) {
		ebin = j;
		break;
	      }
	  }
	    if (ebin == 0)
	      en = FLUX_EN[id - 1][ebin];
	    else {
	      en = ((FLUX_INT[id - 1][ebin - 1]) * FLUX_EN[id - 1][ebin]
		    - (FLUX_INT[id - 1][ebin]) * FLUX_EN[id - 1][ebin - 1]
		    + (prb) * (FLUX_EN[id - 1][ebin - 1] -
			       FLUX_EN[id - 1][ebin]))
		/ ((FLUX_INT[id - 1][ebin - 1]) - (FLUX_INT[id - 1][ebin]));
	    }
	    
	    
	  }
	//MODE VEM = CENTRAL AND VERTICAL 
	if (theConfig ()->PartMode == "VEM") {
	  azim = 2. * PI * Rand ();
	  theta = 0;
	  xx = 0;
	  yy = 0;
	  zz = STATION_HEIGHT * (1. - eps);
	}
	//MODE SCINTILLATOR ( ACCORDING TOPOSITIONS READ IN CALIB.H FILE)
	else if (theConfig ()->PartMode == "SCINTILLATOR") {
	  Double_t xhigh =
	    (Rand () - 0.5) * XSIZE_SCINT_HIGH + XPOS_SCINT_HIGH;
	  Double_t yhigh =
	    (Rand () - 0.5) * YSIZE_SCINT_HIGH + YPOS_SCINT_HIGH;
	  Double_t xlow = (Rand () - 0.5) * XSIZE_SCINT_LOW + XPOS_SCINT_LOW;
	  Double_t ylow = (Rand () - 0.5) * YSIZE_SCINT_LOW + YPOS_SCINT_LOW;
	  Double_t rad =
	    sqrt ((xhigh - xlow) * (xhigh - xlow) +
		  (yhigh - ylow) * (yhigh - ylow));
	  if (ZPOS_SCINT_HIGH != ZPOS_SCINT_LOW)
	    theta = atan (rad / (ZPOS_SCINT_HIGH - ZPOS_SCINT_LOW));
	  else
	    theta = PI / 2.;

	  if (i < 100)
	    thetamoy += theta * RAD2DEG;
	  if (i == 100)
	    if (theta * RAD2DEG < 2)
	      cerr << " VERTICAL PARTICLES " << thetamoy / 100. << endl;
	    else
	      cerr << " THETA OF PARTICLES = " << (thetamoy /
						   100.) << "  degree " <<
		endl;

	  if (xhigh != xlow)
	    azim = atan ((ylow - yhigh) / (xlow - xhigh));
	  else
	    azim = PI / 2;
	  if (xlow < xhigh)
	    azim += PI;
	  if (azim < 0)
	    azim += 
2 * PI;

 	  if (cos (theta) != 0) {
	    if ((tan (theta) * (ZPOS_SCINT_HIGH - STATION_HEIGHT)) >
		(sqrt(xhigh*xhigh+yhigh*yhigh) - STATION_RADIUS)) {
	      xx =
		xhigh + tan (theta) * (ZPOS_SCINT_HIGH -
				       STATION_HEIGHT) * cos (azim);
	      yy =
		yhigh + tan (theta) * (ZPOS_SCINT_HIGH -
				       STATION_HEIGHT) * sin (azim);
	      zz = STATION_HEIGHT * (1. - eps);
	      // cout<<"ici "<<xx<<" "<<yy<<" "<<azim<<endl;
	      //  cout<<theta*RAD2DEG<<" "<<azim*RAD2DEG<<" "<<xhigh<<" "<<yhigh<<" "<<xlow<<" "<<ylow<<" "<<xx<<" "<<yy<<" "<<zz<<endl;	 
	      ntop++;
	    } else {
	      xx = -STATION_RADIUS * cos (azim);
	      yy = -STATION_RADIUS * sin (azim);
	      zz =  ZPOS_SCINT_HIGH - ((sqrt(xhigh*xhigh+yhigh*yhigh)-
				   STATION_RADIUS)) / tan (theta);
	      // cout<<theta*RAD2DEG<<" "<<azim*RAD2DEG<<" "<<xhigh<<" "<<yhigh<<" "<<xlow<<" "<<ylow<<" "<<xx<<" "<<yy<<" "<<zz<<endl;
	     
	      nside++;
	    }
	  } else {
	    xx = -STATION_RADIUS * cos (azim);
	    yy = -STATION_RADIUS * sin (azim);
	    zz = ZPOS_SCINT_HIGH;
	    nside++;
	    cout<<"big bug"<<endl;
	  }
	  // cout<<xlow<<" "<<xhigh<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<tan(theta)*(ZPOS_SCINT_HIGH-STATION_HEIGHT)<<endl;

	} else if (theConfig ()->PartMode == "HORIZONTAL") {

	  azim = 2. * PI * Rand ();
	  theta = 90. * DEG2RAD;
	  y = 2 * Rand () * STATION_RADIUS - STATION_RADIUS;
	  x =
	    sqrt (STATION_RADIUS * STATION_RADIUS - y * y) * (1. - 5. * eps);
	  //  zz = Rand() * STATION_HEIGHT;
	  zz = STATION_HEIGHT / 2.;
	  xx = -x * cos (azim) + y * sin (azim);
	  yy = -x * sin (azim) - y * cos (azim);

	} else {

	  if (theConfig ()->PartMode == "FIXEDTHETA"){
	    theta = theConfig ()->PartTheta * DEG2RAD;
	    azim = 180* DEG2RAD;
}
	  else if (theConfig ()->PartMode == "RANDOM"){	//random
	    theta = acos (pow (1 - Rand (), 1 / 4.));
	    azim = 2. * PI * Rand ();
	  }
	  
	  d_perp = STATION_RADIUS * cos(theta) + STATION_HEIGHT * sin(theta);
	  d_rand = d_perp * Rand ();
	  //d_rand=STATION_HEIGHT * sin (theta);
	  //on the side 
	  if (d_rand < STATION_HEIGHT * sin (theta)) {
	    y = 2 * Rand () * STATION_RADIUS - STATION_RADIUS;
	    x =
	      sqrt (STATION_RADIUS * STATION_RADIUS - y * y) * (1. -
								5. * eps);
	    zz = Rand () * STATION_HEIGHT;
	    xx = -x * cos (azim) + y * sin (azim);
	    yy = -x * sin (azim) - y * cos (azim);
	    nside++;
	  }
	  // on the top
	  else {
	    double r = sqrt (Rand ()) * STATION_RADIUS;
	    double az = Rand () * 2 * PI;
	    x = r * cos (az);
	    y = r * sin (az);
	    xx = x;
	    yy = y;
	    zz = STATION_HEIGHT * (1. - eps);
	    ntop++;
	  }

	  //    zz = 1.2;
	  //      xx =- STATION_RADIUS * cos(azim) ;
//          yy = -STATION_RADIUS * sin(azim);

	}
	///end of all possible cases
	cx = sin (theta) * cos (azim);
	cy = sin (theta) * sin (azim);
	uxy2 = cx * cx + cy * cy;
	if (uxy2 < .999999)
	  cz = -sqrt (1 - uxy2);
	else
	  cz = -0.001;


	//* Creation of an event of one particle in  tank 0 *//
	if (theConfig ()->PartAddMode == "DOUBLE" && ipart > 0)
	  time = 3000 * Rand ();
	else
	  time = 0;

	part =
	  new ParticleInTank (0, id, id, xx, yy, zz, cx, cy, cz, time, en, w);

	part->DoDetSim (theConfig ()->SimMode);

	if (ipart == 0) {
	  hst = new HitStation ();
	  //gives a time flag to the particle according to the mu distribution
	  double randomnumber = Rand ();
	  //      while(randomnumber<0.4999)
	  while (randomnumber == 0)
	    randomnumber = Rand ();
	  //time0+= - 375000 * log(randomnumber); //mu only
	  time0 += -256000 * log (randomnumber);	//all particles
	  hst->SetTime0 (time0);

	  if(theConfig ()->PartMode == "FIXEDTHETA")
	    event = new Event (id, en, theta, azim,xx,yy,zz);
	  else 
	    event = new Event (id, en, theta, azim);

	  if ((i + 1) % 1000 == 0) {
	    cerr << " particles processed  " << i + 1 << endl;
	    //for(int j=0;j<25;j++) cerr<<"\b";
	  }
	}
	hst->CountPart (id, w, time);

	if (theConfig ()->SimMode == "DETAILED")
	  hst->ReadParticleInTank (part);
	//  else  //fast simulation
	//        {
	//          Int_t nbins = part->fBinContent.size();
	//          for(Int_t nb = 0 ; nb<nbins ; nb++){
	//            Double_t bincontent =  part->fBinContent[nb];
	//            hst->WriteADCBin(id,nb,bincontent);
	//          }
	//        }
	delete part;

      }
      if (theConfig ()->SimMode == "DETAILED")
	hst->DoElecSim ();

      event->AddTank (hst);
      event->CountTanks ();
      theConfig ()->OutBranch->SetAddress (&event);
      theConfig ()->OutTree->Fill ();
      delete event;
      delete hst;


    }				//end of loop on events

    if (theConfig ()->PartMode == "FIXEDTHETA"
	|| theConfig ()->PartMode == "RANDOM"
	|| theConfig ()->PartMode == "SCINTILLATOR") {
      cerr << endl << endl << "Entries in the tank :";
      cerr << " nside = " << nside << " ntop = " << ntop << endl << endl;
    }

  }


  if (theConfig ()->Acceptance) {
    ofstream out ("accept.dat", ios::out);
    out << theConfig ()->NEvents << "  " << nt4 << "  " << nt3 << endl;
    out.close ();
  }
  // cleanings
  theConfig ()->OutRootFile->Write ();
  cerr << endl;
  for (int ii = 0; ii < 80; ii++)
    cerr << "-";
  cerr << endl;
  cerr << endl << "\t  ROOT FILE WRITTEN =  " << theConfig ()->
    OutRootFileName << endl;
  cerr << endl;
  for (int ii = 0; ii < 80; ii++)
    cerr << "-";
  cerr << endl;
  theConfig ()->OutRootFile->Close ();


  return 1;
}
