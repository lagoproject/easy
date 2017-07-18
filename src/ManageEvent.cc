
/*
  ManageEvent.cc

  implementation file for ManageEvent.h

*/

#include "ManageEvent.h"
#include "ShowerParam.h"
#include "Array.h"
#include "EasySim.h"
#include "Calib.h"
#include "TRandom3.h"
#include "Utils.h"
#include <stdlib.h>

extern Array *gArray;

Event *
GenerateEvent (UInt_t evid)
{
  Event *event;
  event = new Event (gShParamP->Primary, gShParamP->Energy,
		     gShParamP->Theta, gShParamP->Phi,
		     theConfig ()->EasAvrg, theConfig ()->EasSpread,
		     theConfig ()->NorAvrg, theConfig ()->NorSpread);
  event->SetDepth (gShParamP->XMax, gShParamP->X0);

  cout << endl << "\tEVENT GENERATION ==>" << endl << endl;
  cout << endl << "\tEvent " << evid << " generated  with core at x = "
    << event->fXCore << " y = " << event->fYCore << endl;
  cout << "\tthat is northing = " << event->fNorCore << " easting = " <<
    event->fEasCore << endl;

  return event;
}

void
DoSampling (Double_t xcore, Double_t ycore, vector < SampStation > *tank_list)
{
#ifndef CALIBONLY
  // just to simplify notation, without making copies
  vector < SampStation > &thistanklist = *(tank_list);

  vector < Int_t > nph, nel, nmu;
  vector < Double_t > eph, eel, emu;
  vector < Int_t > nph_samp, nel_samp, nmu_samp;

  Double_t area_top = PI * STATION_RADIUS * STATION_RADIUS;
  Double_t area_side = 2 * STATION_HEIGHT * STATION_RADIUS;

  Double_t cosphi_sh = gShParamP->CosPhi;
  Double_t sinphi_sh = gShParamP->SinPhi;
  Double_t costh_sh = gShParamP->CosTheta;
  Double_t sinth_sh = gShParamP->SinTheta;

  thistanklist.clear ();
  nph.clear ();
  nel.clear ();
  nmu.clear ();
  nph_samp.clear ();
  nel_samp.clear ();
  nmu_samp.clear ();
  eph.clear ();
  eel.clear ();
  emu.clear ();

  cout << endl;
  for (int ii = 0; ii < 80; ii++)
    cout << "-";
  cout << endl;
  cout << "\tSampling in tanks started " << endl;

  //* Here we do a pre-selection of tanks
  //* Loop on all tanks of the array
  Int_t nombtank = gArray->Size ();
  for (Int_t ist = 0; ist < nombtank; ist++) {
    Station *st = &gArray->fStationList[ist];
    Double_t xstat = st->fEasting - xcore;
    Double_t ystat = st->fNorthing - ycore;
    Double_t zstat = st->fAltitude - 1400;


    Double_t xstat_sf =
      (xstat * cosphi_sh + ystat * sinphi_sh) * costh_sh - zstat * sinth_sh;
    Double_t ystat_sf = -xstat * sinphi_sh + ystat * cosphi_sh;
    //      Double_t zstat_sf = (xstat * cosphi_sh + ystat * sinphi_sh) * sinth_sh+ zstat * costh_sh;
    Double_t rstat_sf = sqrt (xstat_sf * xstat_sf + ystat_sf * ystat_sf);

    if (rstat_sf <= theConfig ()->RMax && rstat_sf >= theConfig ()->RMin) {
      //          st->WriteShowerFrame(xstat_sf,ystat_sf,zstat_sf);
      SampStation sst (st->fId, st->fNorthing, st->fEasting, st->fAltitude);
      thistanklist.push_back (sst);
    }
  }				//end of loop on all tanks

  //* Definition of a new list of pre-selected tanks
  Int_t nombtanksel = thistanklist.size ();
  cout << "\tPreselected array of " << nombtanksel << " tanks" << endl;

  //* Initialization of counters
  for (Int_t ist = 0; ist < nombtanksel; ist++) {
    nph_samp.push_back (0);
    nel_samp.push_back (0);
    nmu_samp.push_back (0);
    nph.push_back (0);
    nel.push_back (0);
    nmu.push_back (0);
    eph.push_back (0);
    eel.push_back (0);
    emu.push_back (0);
  }

  //* Starting the loop on particles of the sampling


  cout << "\t --->>>> Sampling " << gShParamP->
    NParticles << " particles ... ";
  fflush (stdout);
  int pc = 0;
  for (Int_t ip = 0; ip < gShParamP->NParticles; ip++) {
    if ((ip % (gShParamP->NParticles / 100)) == 0) {
      if (!ip) {
	cout.width (3);
	cout << pc++ << "%";
      } else {
	cout << "\b\b\b\b";
	cout.width (3);
	cout << pc++ << "%";
	fflush (stdout);
      }
    }
    PartGrnd_ROOT *part = theConfig ()->InShower->showerGrnd ()->GetPart (ip);

    Int_t idpart = part->Id;

    Double_t x = part->X;
    Double_t y = part->Y;

    Double_t ux = part->UX;
    Double_t uy = part->UY;
    Double_t uz = part->UZ;

    Double_t epart = part->E;
    Double_t wpart = part->W;
    Double_t timepart = part->T;

    Double_t norm = sqrt (ux * ux + uy * uy + uz * uz);
    Double_t cxpart = ux / norm;
    Double_t cypart = uy / norm;
    Double_t czpart = uz / norm;

    Double_t thetapart = acos (czpart);
    Double_t tantheta_part = tan (thetapart);

    Double_t phipart;
    Double_t xpart;
    Double_t ypart;

    if (theConfig ()->PhiRotation != 0) {
      phipart = atan2 (cxpart, cypart) + theConfig ()->PhiRotation * DEG2RAD;
      //new values of cxpart and cypart take into account a possible shower rotation;
      // this is to be tested ( )
      xpart =
	x * cos (theConfig ()->PhiRotation * DEG2RAD) -
	y * sin (theConfig ()->PhiRotation * DEG2RAD);
      ypart =
	x * sin (theConfig ()->PhiRotation * DEG2RAD) +
	y * cos (theConfig ()->PhiRotation * DEG2RAD);
      cxpart = sin (thetapart) * cos (phipart);
      cypart = sin (thetapart) * sin (phipart);
    } else {
      xpart = x;
      ypart = y;
    }

    //coordinates of  particles in the "shower frame"
    Float_t xpart_sf =
      (Float_t) (xpart * cosphi_sh + ypart * sinphi_sh) * costh_sh;
    Float_t ypart_sf = (Float_t) (ypart * cosphi_sh - xpart * sinphi_sh);
    Float_t zpart_sf =
      (Float_t) (xpart * cosphi_sh + ypart * sinphi_sh) * sinth_sh;
    Float_t azipart_sf = (Float_t) atan2 (ypart_sf, xpart_sf);
    Float_t rpart_sf =
      (Float_t) sqrt (xpart_sf * xpart_sf + ypart_sf * ypart_sf);
    Float_t time_plane = (Float_t) (-zpart_sf / CLIGHT);

    //  sampling particles within rmax from shower axis
    if (rpart_sf <= theConfig ()->RMax * (1. + theConfig ()->DeltaR)) {
      for (Int_t ist = 0; ist < nombtanksel; ist++) {
	SampStation *sst = &thistanklist[ist];
	Double_t azistat_sf = thistanklist[ist].Azim_sf ();
	Double_t rstat_sf = thistanklist[ist].R_sf ();
	Double_t dazi = azistat_sf - azipart_sf;
	if (dazi > PI)
	  dazi += -2 * PI;
	if (dazi < -PI)
	  dazi += 2 * PI;

	Double_t dr = theConfig ()->DeltaR * rstat_sf;

	if (abs (rstat_sf - rpart_sf) <= dr
	    && abs (dazi) <= theConfig ()->DeltaAzim) {
	  // counting particles in the local sample of the station                
	  if (abs (idpart) == 1)
	    nph_samp[ist]++;
	  else if (abs (idpart) == 2)
	    nel_samp[ist]++;
	  else if (abs (idpart) == 3)
	    nmu_samp[ist]++;

	  //  sampling area for this station (projected on ground)
	  Double_t area_samp =
	    4 * rstat_sf * dr * theConfig ()->DeltaAzim / costh_sh;

	  //  Poisson law for number of top entries
	  Double_t avrg_n = (wpart * area_top / area_samp);

	  Int_t n_entry_top = (Int_t) Rand ("POISSON", avrg_n);
	  // Poisson law for number of side entries (side area is corrected forincidence angle)
	  avrg_n = (wpart * area_side * (-tantheta_part) / area_samp);

	  Int_t n_entry_side = (Int_t) Rand ("POISSON", avrg_n);
	  Int_t n_entry = n_entry_top + n_entry_side;
	  if (n_entry > 0) {
	    Particle *p =
	      new Particle (idpart, n_entry_top, n_entry_side, time_plane,
			    timepart, epart, cxpart, cypart, czpart);
	    sst->AddEntry (p);
	    delete p;
	  }
	}			//end of part in sampling zone of a tank
      }				//end of loop on tanks
    }				//end of first sampling zone
  }				//end of loops on particle

  cout << endl << endl << "\tEnd of sampling" << endl;
  for (int ii = 0; ii < 80; ii++)
    cout << "-";
  cout << endl;
#endif
}

void
AddRandomParticles (vector < SampStation > *samp_stat)
{
  // just to simplify notation, without making copies
  vector < SampStation > &sst = *(samp_stat);

  unsigned int nsst = sst.size ();
  int nrdtot = 0, nrd;
  for (unsigned int ist = 0; ist < nsst; ist++) {
    nrd = (int) Rand ("POISSON", 0.0576);	// muon random rate - 3 kHz
    nrdtot += nrd;
    for (int j = 0; j < nrd; j++) {
      double theta = acos (pow (1. - Rand (), 1. / 4.));
      double azim = 2. * M_PI * Rand ();
      double cx = sin (theta) * cos (azim);
      double cy = sin (theta) * sin (azim);
      double cz;
      if (cx * cx + cy * cy < .999999)
	cz = -sqrt (1. - cx * cx - cy * cy);
      else
	cz = -0.001;
      double time = -6400.+ 19200. * Rand ();
      double d_perp =
	STATION_RADIUS * cos (theta) + STATION_HEIGHT * sin (theta);
      double d_rd = d_perp * Rand ();
      int nside = 0, ntop = 0;
      if (d_rd < STATION_HEIGHT * sin (theta)) {
	nside = 1;
      } else {
	ntop = 1;
      }
      int ebin = 0;
      double prob = Rand ();
      for (int k = 0; k < FLUX_SIZE; k++) {
	if (prob < FLUX_INT[2][k]) {
	  ebin = k;
	  break;
	}
      }
      double ener;
      if (ebin == 0)
	ener = FLUX_EN[2][ebin];
      else {
	ener =
	  (FLUX_INT[2][ebin - 1] * FLUX_EN[2][ebin] -
	   FLUX_INT[2][ebin] * FLUX_EN[2][ebin - 1] +
	   prob * (FLUX_EN[2][ebin - 1] -
		   FLUX_EN[2][ebin])) / (FLUX_INT[2][ebin - 1] -
					 FLUX_INT[2][ebin]);
      }
      Particle *p = new Particle (3, ntop, nside, 0., time, ener, cx, cy, cz);
      sst[ist].AddEntry (p);
      delete p;
    }
  }


  cout << "\tAdded " << nrdtot << " random particles." << endl << endl;
  for (int ii = 0; ii < 80; ii++)
    cout << "-";
  cout << endl;
}


HitStation *
GenerateHitStation (const SampStation & sst)
{
  int nombpart = sst.fPartList.size ();

  Double_t cosphi_sh = gShParamP->CosPhi;
  Double_t sinphi_sh = gShParamP->SinPhi;
  Double_t costh_sh = gShParamP->CosTheta;
  Double_t sinth_sh = gShParamP->SinTheta;

  HitStation *hst =
    new HitStation (sst.fId, sst.fNorthing, sst.fEasting, sst.fAltitude);

  if (nombpart <= 0) {
    return hst;
  }

  cout << "\tSimulation of tank " << hst->fId << " ... ";
  hst->DoDetSim (theConfig ()->SimMode, costh_sh, sinth_sh, cosphi_sh,
		 sinphi_sh, sst);
  cout << endl;

  int t2trig = 0;
  if (theConfig ()->SimMode == "DETAILED") {
    hst->DoElecSim ();
    t2trig = hst->DoTrigSim ();
    hst->SetTime0 ();
    cout << "\t\tDistance from the core in the shower frame : " <<
      floor (hst->R_sf ()) << " m" << endl;
    cout << "\t\tSignal : " << floor (hst->SignalInVem () * 100) /
      100. << " VEM" << endl;
    cout << "\t\tParticles in tank --- total number :" << hst->
      fNPartTot << endl;
    cout << "\t\t                  --- #muons       :" << hst->fNmu << endl;
    cout << "\t\t                  --- #photons     :" << hst->fNph << endl;
    cout << "\t\t                  --- #e+/e-       :" << hst->fNel << endl;
    cout << "\t\tPhoto-electrons   --- total number :" << floor (hst->
								 fNpe) <<
      endl;
    if (theConfig ()->MuEmMode == "MUEM") {
      cout << "\t\t                  --- #pe muons    :" << floor (hst->
								   fNpe_mu) <<
	endl;
      cout << "\t\t                  --- #pe photons  :" << floor (hst->
								   fNpe_ph) <<
	endl;
      cout << "\t\t                  --- #pe e+/e-    :" << floor (hst->
								   fNpe_el) <<
	endl;
    }
  }


  if (t2trig) {
    cout << "\t\t=> Station triggered the T2 level" << endl;
  } else if (hst->fT1Threshold) {
    cout << "\t\t=> Station triggered the T1 level" << endl;
  } else {
    cout << "\t\t=> Station didn't pass any trigger" << endl;
  }

  return hst;
}
