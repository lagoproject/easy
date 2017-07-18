#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <sys/types.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>
// #include "ShowerAnalyze.h"
#include "Event.h"
#include "Trigger.h"
#include "ElecConstants.h"


#include "IoSd.h"
#include "Ec.h"
TROOT root ("Main", "PAO Simulation to CDAS");	// ROOT workspace

int
main (int argc, char *argv[])
{
  TFile *inrootfile;
  char outrootfilename[256];
  char *tmp;
   for (int ifile = 1; ifile < argc; ifile++) {

    inrootfile = new TFile (argv[ifile], "READ");	// Open file
    tmp = strstr (argv[ifile], "_");
    strcpy (outrootfilename, "cdas");
    strcat (outrootfilename, tmp);
    cout << outrootfilename << endl;
    IoSd outiosd (outrootfilename, "w");

    Event *event = new Event;
    TTree *intree = (TTree *) inrootfile->Get ("TreeEvent");	// Get Tree
    intree->Print ();
    TBranch *inbranch = intree->GetBranch ("event");
    inbranch->SetAddress (&event);

    int nbentries = (int) intree->GetEntries ();
    int *flagtrig;

    for (int itirage = 0; itirage < nbentries; itirage++) {
      IoSdEvent cdasev;
      intree->GetEvent (itirage);
      kEc::StrictCalib(0);
      cdasev.Id = itirage;
      cdasev.Simul = new IoSdSimul ();
      cdasev.Simul->Primary = event->fPrimary;
      cdasev.Simul->Theta = event->fTheta;
      cdasev.Simul->Phi = event->fAzim;
      cdasev.Simul->E = event->fEnergy;
      cdasev.Simul->XMax = event->fXmax;
      cdasev.Simul->XFirst = event->fX0;
      cdasev.Simul->CoreEasting = event->fEasCore;
      cdasev.Simul->CoreNorthing = event->fNorCore;
      strcpy (cdasev.Simul->MotherProgram, "Corsika v??");
      strcpy (cdasev.Simul->ShowerName, "???");
      cdasev.Trigger.Id = 0;
      cdasev.Trigger.PreviousId = 0;
      cdasev.Trigger.NumberOfStation = event->fNombTank;
      cdasev.Trigger.Second = 1072961605 - (315964800 - 13);
      cdasev.Trigger.MicroSecond = 0;
      cdasev.Trigger.Sender = "EASYSIM";
      if (event->fT4Algo == 1)
	cdasev.Trigger.Algo = "3TOT";
      else if (event->fT4Algo == 2)
	cdasev.Trigger.Algo = "4C1";
      else
	cdasev.Trigger.Algo = "NADA";
      cdasev.Trigger.SDPAngle = 0.;

      vector < HitStation > *hstlist = &event->fHitStationList;
      int nstattrig = 0;
      flagtrig = new int[event->fNombTank];
      for (int is = 0; is < event->fNombTank; is++) {
	if ((*hstlist)[is].fT2Threshold || (*hstlist)[is].fT2ToT || (*hstlist)[is].fT1Threshold) {
	  nstattrig++;
	  flagtrig[is] = 0;
	} else
	  flagtrig[is] = 1;
      }

      cdasev.NumberOfErrorZeroStation = nstattrig;
      cdasev.Stations.resize (event->fNombTank);

      for (int is = 0; is < event->fNombTank; is++) {
	if (!(*hstlist)[is].fT2Threshold && !(*hstlist)[is].fT2ToT && !(*hstlist)[is].fT1Threshold) {
	  cdasev.Stations[is].Error = 2;
	  cdasev.Stations[is].Id = (*hstlist)[is].fId;
	  cdasev.Stations[is].Name = "";
	  cdasev.Stations[is].Easting = (*hstlist)[is].fEasting;
	  cdasev.Stations[is].Northing = (*hstlist)[is].fNorthing;
	  cdasev.Stations[is].Altitude = (*hstlist)[is].fAltitude;
	  cdasev.Stations[is].Trigger.Window = 30;
	  continue;
	}
	cdasev.Stations[is].Error = 0;
	cdasev.Stations[is].Id = (*hstlist)[is].fId;
	cdasev.Stations[is].Name = "";
	cdasev.Stations[is].Easting = (*hstlist)[is].fEasting;
	cdasev.Stations[is].Northing = (*hstlist)[is].fNorthing;
	cdasev.Stations[is].Altitude = (*hstlist)[is].fAltitude;
	if ((*hstlist)[is].fT2ToT) {
	  cdasev.Stations[is].Trigger.Type = 2;
	  cdasev.Stations[is].Trigger.Name = "TOT";
	  cdasev.Stations[is].Trigger.Window = 0;
	} else if ((*hstlist)[is].fT2Threshold) {
	  cdasev.Stations[is].Trigger.Type = 1;
	  cdasev.Stations[is].Trigger.Name = "Threshold";
	  cdasev.Stations[is].Trigger.Window = 0;
	} else if ((*hstlist)[is].fT1Threshold) {
	  cdasev.Stations[is].Trigger.Type = 1;
	  cdasev.Stations[is].Trigger.Name = "Threshold";
	  cdasev.Stations[is].Trigger.Window = 30;
	}

	cdasev.Stations[is].Trigger.Offset = 0;

	cdasev.Stations[is].Gps = new IoSdGps ();
	cdasev.Stations[is].Gps->Second = 1072961605 - (315964800 - 13);
	cdasev.Stations[is].Gps->Tick =
	  int ((500000000. + (*hstlist)[is].fTime0) / 10.);
	cdasev.Stations[is].Gps->Current100 =
	  cdasev.Stations[is].Gps->Next100 = 100000000;
	cdasev.Stations[is].Gps->Current40 = cdasev.Stations[is].Gps->Next40 =
	  40000000;
	cdasev.Stations[is].Gps->PreviousST =
	  cdasev.Stations[is].Gps->CurrentST = 0;
	cdasev.Stations[is].Gps->NextST = cdasev.Stations[is].Gps->Offset = 0;
	// bof...
	cdasev.Stations[is].Gps->CurrentST =
	  (((int) (*hstlist)[is].fTime0 % 10) == 0) ? 0 : 5;
	//      cout<<cdasev.Stations[is].Gps->Tick<<"  "<<(*hstlist)[is].fTime0<<"  "<<cdasev.Stations[is].Gps->CurrentST<<"  "<<(int) ((cdasev.Stations[is].Gps->Tick*(1000000000.0/cdasev.Stations[is].Gps->Next100)+cdasev.Stations[is].Gps->CurrentST))<<endl;

	cdasev.Stations[is].Fadc = new IoSdFadc ();
	cdasev.Stations[is].Fadc->NSample = 768;
	for (int ipm = 0; ipm < NPM; ipm++) {
	  for (int ig = 0; ig < 2; ig++) {
	    for (int ibin = 0; ibin < 768; ibin++) {
	      cdasev.Stations[is].Fadc->Trace[ipm][ig][ibin] = (short) 50;
	    }
	    int deltabin = ((*hstlist)[is].fTriggerBin >= 255) ? 255 : (*hstlist)[is].fTriggerBin;
	    for (int ibin = -deltabin; ibin < MAXNUMBEROFADCBINS - 255; ibin++) {
	      cdasev.Stations[is].Fadc->Trace[ipm][ig][255 + ibin] +=
		(*hstlist)[is].fPMT[ipm].fADC[1 - ig][ibin + (*hstlist)[is].fTriggerBin];
	    }
	  }
	}

	/*
	   for (int  ipm=0; ipm<3; ipm++) {
	   for (int ig=0; ig<2; ig++) {
	   for (int ibin=0; ibin<256; ibin++) {
	   cdasev.Stations[is].Fadc->Trace[ipm][ig][ibin] =  (short) 0;
	   }                
	   for (int ibin=0; ibin<MAXNUMBEROFADCBINS-firstADCbin; ibin++) {
	   cdasev.Stations[is].Fadc->Trace[ipm][ig][ibin+256] =  
	   (*hstlist)[is].fPMT[ipm].fADC[1-ig][firstADCbin+ibin];             
	   //cout<<(*hstlist)[is].fPMT[ipm].fADC[1-ig][firstADCbin+ibin]<<"  ";
	   //cdasev.Stations[is].Fadc->Trace[ipm][ig][ibin] =  (*hstlist)[is].fPMT[ipm].fADC[1-ig][ibin];
	   }
	   //       cout<<endl;
	   for (int ibin=MAXNUMBEROFADCBINS-firstADCbin+256; ibin<768; ibin++) {
	   cdasev.Stations[is].Fadc->Trace[ipm][ig][ibin] = 0;
	   }
	   }
	   }
	 */

	cdasev.Stations[is].Calib = new IoSdCalib ();
	cdasev.Stations[is].Calib->Version = 4;
	cdasev.Stations[is].Calib->TubeMask = 7;
	cdasev.Stations[is].Calib->StartSecond = 1;
	cdasev.Stations[is].Calib->EndSecond = 3;
	cdasev.Stations[is].Calib->NbT1 = 100;
	cdasev.Stations[is].Calib->NbT2 = 20;
	cdasev.Stations[is].Calib->NbTOT = 1;

	for (int ipm = 0; ipm < NPM; ipm++) {
	  cdasev.Stations[is].Calib->TubeOk[ipm] = 1;
	  cdasev.Stations[is].Calib->Base[ipm] = 50;
	  cdasev.Stations[is].Calib->SigmaBase[ipm] = 0;
	  cdasev.Stations[is].Calib->Base[ipm + NPM] = 50;
	  cdasev.Stations[is].Calib->SigmaBase[ipm + NPM] = 0;
	  cdasev.Stations[is].Calib->NbTDA[ipm] = 0;
	  cdasev.Stations[is].Calib->Evolution[ipm] = 0;
	  cdasev.Stations[is].Calib->VemCharge[ipm] = VEMCHARGEVALUEINADC;
	  cdasev.Stations[is].Calib->VemPeak[ipm] = VEMPEAKVALUEINADC;
	  cdasev.Stations[is].Calib->VemCharge[ipm] /= 1.04;
	  cdasev.Stations[is].Calib->VemPeak[ipm] /= 1.01;
	  // balance correction from Laura
	  cdasev.Stations[is].Calib->VemCharge[ipm] *= 1.01;
	  cdasev.Stations[is].Calib->VemPeak[ipm] *= 0.87;
	  cdasev.Stations[is].Calib->Rate[ipm] = 7000;
	  cdasev.Stations[is].Calib->DA[ipm] = DYNODETOANODE;
	  //Calib->SigmaDA[i]=  DYNODEAMP*DYNODEAMP ;
	  cdasev.Stations[is].Calib->SigmaDA[ipm] = 0;
	}

      }
     
      outiosd.Write (cdasev);

      for (int is = 0; is < event->fNombTank; is++) {
	if (flagtrig[is]) {
	  delete cdasev.Stations[is].Gps;
	  delete cdasev.Stations[is].Calib;
	  delete cdasev.Stations[is].Fadc;
	}
      }
      delete flagtrig;

    }

    delete inrootfile;
    delete event;

    outiosd.Close ();
  }

  return 1;
}
