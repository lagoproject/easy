/*
  Station.cc
  implementation file for all classes Station
*/

#include "Station.h"
#include "Array.h"
#include "Trigger.h"
#include "ElecConstants.h"
#include "ShowerParam.h"
#include "EasySim.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <strstream>
#include "Utils.h"


// ClassImp(Station)
// ClassImp(SampStation)
// ClassImp(PMT)
// ClassImp(HitStation)

//----------------------------------------------------------------
/*
  class Station (parent class)
  A Station is defined by its coordinates.
*/
//----------------------------------------------------------------

Station::Station ()
{

}

Station::Station (Int_t id, Float_t north, Float_t east, Float_t alt)
{
  fId = id;
  fNorthing = north;
  fEasting = east;
  fAltitude = alt;
  fX_sf = X_sf ();
  fY_sf = Y_sf ();
  fZ_sf = Z_sf ();
  fR_sf = R_sf ();
  fAzim_sf = Azim_sf ();
}

Station::Station (Int_t id, Float_t north, Float_t east)
{
  fId = id;
  fNorthing = north;
  fEasting = east;
  fAltitude = ALTITUDE_SITE;
  fX_sf = X_sf ();
  fY_sf = Y_sf ();
  fZ_sf = Z_sf ();
  fR_sf = R_sf ();
  fAzim_sf = Azim_sf ();
}

Station::~Station ()
{

}


// void Station::WriteName(string name)
// {
//   fName=name;
// }

// void Station::WriteShowerFrame(Double_t xstat,Double_t ystat,Double_t zstat)
// {
//   fX_sf=xstat;
//   fY_sf=ystat;
//   fZ_sf=zstat;
//   fR_sf=sqrt(fX_sf * fX_sf + fY_sf * fY_sf );
//   fAzim_sf=atan2(fY_sf,fX_sf);

// }


double
Station::X ()
{
  double x, y;
  UTMToXY (fNorthing, fEasting, &x, &y);
  return x;
}

double
Station::Y ()
{
  double x, y;
  UTMToXY (fNorthing, fEasting, &x, &y);
  return y;
}

double
Station::Z ()
{
  return fAltitude;
}

double
Station::X_sf ()
{
  Double_t cosphi_sh = gShParamP->CosPhi;
  Double_t sinphi_sh = gShParamP->SinPhi;
  Double_t costh_sh = gShParamP->CosTheta;
  Double_t sinth_sh = gShParamP->SinTheta;
  Double_t xstat = fEasting - gShParamP->EasCore;
  Double_t ystat = fNorthing - gShParamP->NorCore;
  Double_t zstat = fAltitude - 1400;

  return (xstat * cosphi_sh + ystat * sinphi_sh) * costh_sh -
    zstat * sinth_sh;
}

double
Station::Y_sf ()
{
  Double_t cosphi_sh = gShParamP->CosPhi;
  Double_t sinphi_sh = gShParamP->SinPhi;
  Double_t xstat = fEasting - gShParamP->EasCore;
  Double_t ystat = fNorthing - gShParamP->NorCore;

  return -xstat * sinphi_sh + ystat * cosphi_sh;
}

double
Station::Z_sf ()
{
  if (theConfig ()->Mode == "CALIB")
    return (0);
  Double_t cosphi_sh = gShParamP->CosPhi;
  Double_t sinphi_sh = gShParamP->SinPhi;
  Double_t costh_sh = gShParamP->CosTheta;
  Double_t sinth_sh = gShParamP->SinTheta;
  Double_t xstat = fEasting - gShParamP->EasCore;
  Double_t ystat = fNorthing - gShParamP->NorCore;
  Double_t zstat = fAltitude - 1400;

  return (xstat * cosphi_sh + ystat * sinphi_sh) * sinth_sh +
    zstat * costh_sh;
}

double
Station::R_sf ()
{
  double xsf = X_sf (), ysf = Y_sf ();	//, zsf=Z_sf();
  return sqrt (xsf * xsf + ysf * ysf);
}

double
Station::Azim_sf ()
{
  double xsf = X_sf (), ysf = Y_sf ();
  return atan2 (ysf, xsf);
}

//----------------------------------------------------------------
/*
  class SampStation
  A SampStation contains all particles that have to be simulated
*/
//----------------------------------------------------------------


SampStation::SampStation ()
{

}

SampStation::SampStation (Int_t id, Float_t north, Float_t east, Float_t alt)
:Station (id, north, east, alt)
{
  fNmu = 0;
  fNph = 0;
  fNel = 0;
  fNPartTot = 0;
}

SampStation::SampStation (Int_t id, Float_t north, Float_t east)
:Station (id, north, east)
{
  fNmu = 0;
  fNph = 0;
  fNel = 0;
  fNPartTot = 0;
}

SampStation::~SampStation ()
{

}

void
SampStation::AddEntry (Particle * p)
{
  fPartList.push_back (*p);
  Int_t w = (Int_t) (p->fWtop + p->fWside);
  if (p->fId == 1 || p->fId == -1)
    fNph += w;
  if (p->fId == 2 || p->fId == -2)
    fNel += w;
  if (p->fId == 3 || p->fId == -3)
    fNmu += w;
  fNPartTot += w;
}

//----------------------------------------------------------------
/*
  class PMT
  A PMT.
*/
//----------------------------------------------------------------

PMT::PMT ()
{
  fNpe = 0;
  fNpe_direct = 0;
  fNpe_mu = 0;
  fNpe_ph = 0;
  fNpe_el = 0;
  fFirstPelTime = 0;
  fFirstPelBin = MAXNUMBEROFTIMEBINS;
  fFirstADCBin = 0;
  fSampFact = 1.;

  for (Int_t i = 0; i < MAXNUMBEROFTIMEBINS; i++) {
    fFEFilter_hi[i] = 0;
    fFEFilter_lo[i] = 0;
  }

  for (Int_t i = 0; i < MAXNUMBEROFADCBINS; i++) {
    fADC_ns[i] = 0;
    for (Int_t j = 0; j < 2; j++) {
      fADC[j][i] = 0;
      fADC_mu[j][i] = 0;
      fADC_em[j][i] = 0;
    }
  }

}

PMT::~PMT ()
{

}

// PMT::PMT(const PMT& p)
// {
//   *this=p;
// }

void
PMT::Init ()
{
  fNpe = 0;
  fNpe_direct = 0;
  fNpe_mu = 0;
  fNpe_ph = 0;
  fNpe_el = 0;
  fFirstPelTime = 0;
  fFirstPelBin = MAXNUMBEROFTIMEBINS;
  fFirstADCBin = 0;
  fSampFact = 1.;
  fTimeProfile.clear ();
  fTimeProfile_mu.clear ();
  fTimeProfile_em.clear ();
  fPMTSignal_hi.clear ();
  fPMTSignal_lo.clear ();


  for (Int_t i = 0; i < MAXNUMBEROFTIMEBINS; i++) {
    fFEFilter_hi[i] = 0;
    fFEFilter_lo[i] = 0;
  }

  for (Int_t i = 0; i < MAXNUMBEROFADCBINS; i++)
    for (Int_t j = 0; j < 2; j++) {
      fADC[j][i] = 0;
      fADC_mu[j][i] = 0;
      fADC_em[j][i] = 0;
    }
}

void
PMT::WritePE (Int_t idpart, Int_t time)
{

  if (fTimeProfile.count (time))
    fTimeProfile[time]++;
  else
    fTimeProfile[time] = 1;
  fNpe += (1. / fSampFact);

  if (theConfig ()->MuEmFlag) {

    if (idpart == 1 || idpart == -1) {

      if (fTimeProfile_em.count (time))
	fTimeProfile_em[time]++;
      else
	fTimeProfile_em[time] = 1;
      fNpe_ph += (1. / fSampFact);


    } else if (idpart == 2 || idpart == -2) {
      if (fTimeProfile_em.count (time))
	fTimeProfile_em[time]++;
      else
	fTimeProfile_em[time] = 1;
      fNpe_el += (1. / fSampFact);

    } else if (idpart == 3 || idpart == -3) {

      if (fTimeProfile_mu.count (time))
	fTimeProfile_mu[time]++;
      else
	fTimeProfile_mu[time] = 1;
      fNpe_mu += (1. / fSampFact);
    }
  }
}

void
PMT::DirectLight (int npe_direct)
{
  fNpe_direct += npe_direct;
}

void
PMT::DoElecSim ()
{
  if (theConfig ()->ElecMode == "FULL") {


    PMTAmplification (fTimeProfile, fPMTSignal_lo, fPMTSignal_hi);

    FEFilter (fPMTSignal_hi, fFEFilter_hi);
    FEFilter (fPMTSignal_lo, fFEFilter_lo);

    FADCSampling (fFEFilter_hi, fADC[0]);
    FADCSampling (fFEFilter_lo, fADC[1]);
  }

  else {
    map < Int_t, Double_t > pmtsignal_hi, pmtsignal_lo;
    Double_t fefilter_hi[MAXNUMBEROFTIMEBINS],
      fefilter_lo[MAXNUMBEROFTIMEBINS];
    for (Int_t i = 0; i < MAXNUMBEROFTIMEBINS; i++) {
      fefilter_hi[i] = 0;
      fefilter_lo[i] = 0;

    }

    PMTAmplification (fTimeProfile, pmtsignal_lo, pmtsignal_hi);

    FEFilter (pmtsignal_hi, fefilter_hi);
    FEFilter (pmtsignal_lo, fefilter_lo);

    if (theConfig ()->ElecMode != "SHOWSAT")
      FADCSampling (fefilter_lo, fADC[1]);
    else
      FADCSamplingSat (fefilter_lo, fADC[1], fADC_ns);

    FADCSampling (fefilter_hi, fADC[0]);	//this order to get the firstbin for high gain channel
    if (theConfig ()->ElecMode == "PM") {
      fPMTSignal_hi = pmtsignal_hi;
      fPMTSignal_lo = pmtsignal_lo;
    }

  }

  if (theConfig ()->MuEmFlag) {
    map < Int_t, Double_t > pmtsignal_mu_hi, pmtsignal_mu_lo;
    map < Int_t, Double_t > pmtsignal_em_hi, pmtsignal_em_lo;

    Double_t fefilter_mu_hi[MAXNUMBEROFTIMEBINS],
      fefilter_mu_lo[MAXNUMBEROFTIMEBINS];
    Double_t fefilter_em_hi[MAXNUMBEROFTIMEBINS],
      fefilter_em_lo[MAXNUMBEROFTIMEBINS];

    for (Int_t i = 0; i < MAXNUMBEROFTIMEBINS; i++) {
      fefilter_mu_hi[i] = 0;
      fefilter_em_hi[i] = 0;
      fefilter_mu_lo[i] = 0;
      fefilter_em_lo[i] = 0;
    }

    PMTAmplification (fTimeProfile_mu, pmtsignal_mu_lo, pmtsignal_mu_hi);
    PMTAmplification (fTimeProfile_em, pmtsignal_em_lo, pmtsignal_em_hi);


    FEFilter (pmtsignal_mu_hi, fefilter_mu_hi);
    FEFilter (pmtsignal_em_hi, fefilter_em_hi);
    FEFilter (pmtsignal_mu_lo, fefilter_mu_lo);
    FEFilter (pmtsignal_em_lo, fefilter_em_lo);
    FADCSampling (fefilter_mu_hi, fADC_mu[0]);
    FADCSampling (fefilter_mu_lo, fADC_mu[1]);
    FADCSampling (fefilter_em_hi, fADC_em[0]);
    FADCSampling (fefilter_em_lo, fADC_em[1]);
  }
}


void
PMT::PMTAmplification (map < Int_t, Double_t > &spec, map < Int_t,
		       Double_t > &specAn, map < Int_t, Double_t > &specDy)
{

  if (spec.size () == 0) {
    //  cerr << " PmAmplification :: empty photoelectrons " << endl;
    return;			// empty input
  }


  map < Int_t, Double_t >::iterator itSlot;
  map < Int_t, Double_t >::iterator start = spec.begin ();
  map < Int_t, Double_t >::iterator stop = spec.end ();
  specAn.clear ();
  specDy.clear ();


  //here we consider some transit time fluctuations

  // loop on each photoelectron TimeSlot

  for (itSlot = start; itSlot != stop; itSlot++) {
    Int_t tBin = itSlot->first;	// time bin on photoelectrons channels
    Int_t nPe = (Int_t) itSlot->second;	// number of  photoelectrons for this bin
    //loop on pe's in this time bin
    for (Int_t iPe = 0; iPe < nPe; iPe++) {
      Int_t chargebin =
	(Int_t) (rint (Rand () * (PHOTOELCHARGEPROBSIZE - 2)));
      Double_t scaling =
	PHOTOELCHARGEPROB[chargebin] +
	Rand () * (PHOTOELCHARGEPROB[chargebin + 1] -
		   PHOTOELCHARGEPROB[chargebin]);
      Int_t tBintr = (Int_t) (rint (Rand ("GAUSS", tBin, PMTRANSITTIME)));

      for (Int_t shBin = 0; shBin < PHOTOELPULSESHAPESIZE; shBin++) {
	Int_t iBin =
	  tBintr + PHOTOELPULSESHAPETIME[shBin] + PHOTOELPULSESTARTBIN;
	if (iBin < MAXNUMBEROFTIMEBINS) {
	  specAn[iBin] += 0.001 * scaling * PMGAIN / PMGAINFORPULSESHAPE *
	    PHOTOELPULSESHAPE[shBin] / DYNODETOANODE;
	  specDy[iBin] += 0.001 * scaling * PMGAIN / PMGAINFORPULSESHAPE *
	    PHOTOELPULSESHAPE[shBin];

	}
      }

    }
  }

}


//------------------------------------------------------------------
/*! \fn void SDElecSim::FEFilter(map<Int_t,Double_t>& spec1,map<Int_t,Double_t>& spec2)
    \brief This function does the Pm current filtering in a kind of RC filter
    \param spec1 is the current time pattern take for input argument
    \param spec2 is the voltage time pattern take for output argument
 */
//------------------------------------------------------------------

void
PMT::FEFilter (map < Int_t, Double_t > &spec1, Double_t * spec2)
{


  Double_t *x;
  Double_t *y;
  map < Int_t, Double_t >::iterator beginmap = spec1.begin ();
  map < Int_t, Double_t >::iterator endmap = spec1.end ();
  Int_t firstbin, lastbin;

  Int_t i, bin;
  // static Double_t mem[5];

  if (!spec1.size ()) {

    return;
  }

  endmap--;
  lastbin = endmap->first;
  firstbin = beginmap->first;
  if (firstbin < 0)
    firstbin = 0;

  x = new Double_t[MAXNUMBEROFTIMEBINS];
  for (bin = 0; bin < MAXNUMBEROFTIMEBINS; bin++)
    if (spec1.count (bin))
      x[bin] = spec1[bin];
    else
      x[bin] = 0;

  y = new Double_t[MAXNUMBEROFTIMEBINS + 5];
  for (i = 0; i < MAXNUMBEROFTIMEBINS + 5; i++) {
    y[i] = 0;
  }

  /* Filter recursive equation */


  for (i = firstbin + 5; i < lastbin + 5; i++) {
    y[i] =
      -FEDCOEF[0] * x[i - 5] + FEDCOEF[1] * y[i - 1] + FEDCOEF[2] * y[i - 2]
      + FEDCOEF[3] * y[i - 3] + FEDCOEF[4] * y[i - 4] + FEDCOEF[5] * y[i - 5];


  }


  for (i = firstbin + 5; i < lastbin + 5; i++) {
    y[i] *= 1 / fSampFact;
  }

  memcpy ((void *) spec2, (void *) (y + 5),
	  MAXNUMBEROFTIMEBINS * sizeof (Double_t));



  delete[]y;
  delete[]x;

}


//---------------------------------------------------------------------
/*! \fn void SDElecSim::FADCSampling(Double_t* spec1 ,map<Int_t,Int_t>& spec2)
    \brief This function does the digitization of a FADC
    \param spec1 is the input voltage time pattern
    \param spec2 is the output digitized time pattern
 */
//---------------------------------------------------------------------

void
PMT::FADCSampling (Double_t * spec1, Short_t * spec2)
{
  Int_t i, j, n;
  Double_t value;
  short firstbinfound = 0;

  /* Initial random phase */
  j = (Int_t) (Rand () * (SAMPLINGTIME - 1));

  for (i = 0; i < MAXNUMBEROFADCBINS /*-TIMEOFFSET/SAMPLINGTIME*/ ; i++) {
    value =
      (*(spec1 + j + (i * SAMPLINGTIME) /*+TIMEOFFSET */ ) *
       (Double_t) NBFADCCHANNEL);
    n = (Int_t) value;
    if (Rand () < value - n)
      n++;
    if (n < 0)
      n = 0;
    if (n > NBFADCCHANNEL)
      n = NBFADCCHANNEL;
    spec2[i] = (Short_t) n;
    if (!firstbinfound && spec2[i] >= 1) {
      fFirstADCBin = i;
      firstbinfound = 1;

    }
  }

}


//---------------------------------------------------------------------
/*! \fn void SDElecSim::FADCSamplingsat(Double_t* spec1 ,map<Int_t,Int_t>& spec2)
    \brief This function does the digitization of a FADC
    \param spec1 is the input voltage time pattern
    \param spec2 is the output digitized time pattern
    \ one extra output for  with no saturation
 */
//---------------------------------------------------------------------

void
PMT::FADCSamplingSat (Double_t * spec1, Short_t * spec2, Int_t * spec3)
{

  Int_t i, j, n;
  Double_t value;
  short firstbinfound = 0;


  /* Initial random phase */
  j = (Int_t) (Rand () * (SAMPLINGTIME - 1));

  for (i = 0; i < MAXNUMBEROFADCBINS; i++) {
    value = *(spec1 + j + (i * SAMPLINGTIME)) * (Double_t) NBFADCCHANNEL;
    n = (Int_t) value;
    if (Rand () < value - n)
      n++;
    if (n < 0)
	    n = 0;
    if (n > NBFADCCHANNEL)
	     n = NBFADCCHANNEL;
    spec2[i] = (Short_t) n;
    if (!firstbinfound && spec2[i] >= 1) {
	    fFirstADCBin = i;
	    firstbinfound = 1;
    }
    spec3[i] = n; 

  }
}



//----------------------------------------------------------------
/*
  class HitStation
  A HitStation contains all informations after simulation
*/
//----------------------------------------------------------------


HitStation::HitStation ()
{
  fAltitude = ALTITUDE_SITE;
  fSampFact = 1.;
  fMuTimes.clear ();
  fTimeProfile.clear ();
  fTimeProfile_mu.clear ();
  fTimeProfile_em.clear ();

  for (Int_t i = 0; i < NPM; i++)
    fPMT[i].Init ();

}

HitStation::HitStation (Int_t id, Float_t north, Float_t east, Float_t alt)
  //  : SampStation(id, north, east, alt)
{
  fId = id;
  fNorthing = north;
  fEasting = east;
  fAltitude = alt;
  fX_sf = X_sf ();
  fY_sf = Y_sf ();
  fZ_sf = Z_sf ();
  fR_sf = R_sf ();
  fAzim_sf = Azim_sf ();
  fFirstPelTime = 0;
  fFirstPelBin = MAXNUMBEROFTIMEBINS;
  fT2Threshold = fT1Threshold = fT2ToT = 0;
  fTime0 = -Z_sf () / CLIGHT;
  fTriggerBin = 0;
  fNpe = fNpe_direct = fNpe_mu = fNpe_ph = fNpe_el = 0;
  fNtop = fNside = 0;
  fSampFact = 1.;
  fMuTimes.clear ();
  fTimeProfile.clear ();
  fTimeProfile_mu.clear ();
  fTimeProfile_em.clear ();

  for (Int_t i = 0; i < MAXNUMBEROFADCBINS; i++) {
    fADC_ns[i] = 0;
    for (Int_t j = 0; j < 2; j++) {
      fADC[j][i] = fADC_mu[j][i] = fADC_em[j][i] = 0;
      fDoubleADC[j][i] = fDoubleADC_mu[j][i] = fDoubleADC_em[j][i] = 0;
    }
  }
  for (Int_t i = 0; i < NPM; i++)
    fPMT[i].Init ();

}

HitStation::HitStation (Int_t id, Float_t north, Float_t east)
:SampStation (id, north, east)
{
  fId = id;
  fNorthing = north;
  fEasting = east;
  fX_sf = X_sf ();
  fY_sf = Y_sf ();
  fZ_sf = Z_sf ();
  fR_sf = R_sf ();
  fAzim_sf = Azim_sf ();
  fNbOfBinsOverThres = 0;
  fFirstPelTime = 0;
  fFirstPelBin = MAXNUMBEROFTIMEBINS;
  fT2Threshold = fT1Threshold = fT2ToT = 0;
  fTime0 = -Z_sf () / CLIGHT;
  fTriggerBin = 0;
  fNpe = fNpe_direct = fNpe_mu = fNpe_ph = fNpe_el = 0;
  fNtop = fNside = 0;
  fSampFact = 1.;
  fMuTimes.clear ();
  fTimeProfile.clear ();
  fTimeProfile_mu.clear ();
  fTimeProfile_em.clear ();

  for (Int_t i = 0; i < MAXNUMBEROFADCBINS; i++) {
    fADC_ns[i] = 0;
    for (Int_t j = 0; j < 2; j++) {
      fADC[j][i] = fADC_mu[j][i] = fADC_em[j][i] = 0;
      fDoubleADC[j][i] = fDoubleADC_mu[j][i] = fDoubleADC_em[j][i] = 0;
    }
  }
  for (Int_t i = 0; i < NPM; i++)
    fPMT[i].Init ();

}

HitStation::~HitStation ()
{

}

void
HitStation::WriteSamplingFactor (Double_t sampfact)
{
  fSampFact = sampfact;
  for (Int_t ipm = 0; ipm < NPM; ipm++)
    fPMT[ipm].fSampFact = sampfact;
}

void
HitStation::WritePE (Int_t idpart, Int_t time)
{

  if (fTimeProfile.count (time))
    fTimeProfile[time]++;
  else
    fTimeProfile[time] = 1;
  fNpe += (1. / fSampFact);

  if (theConfig ()->MuEmFlag) {
    if (idpart == 1 || idpart == -1) {
      if (fTimeProfile_em.count (time))
	fTimeProfile_em[time]++;
      else
	fTimeProfile_em[time] = 1;
      fNpe_ph += (1. / fSampFact);

    } else if (idpart == 2 || idpart == -2) {
      if (fTimeProfile_em.count (time))
	fTimeProfile_em[time]++;
      else
	fTimeProfile_em[time] = 1;
      fNpe_el += (1. / fSampFact);

    } else if (idpart == 3 || idpart == -3) {
      if (fTimeProfile_mu.count (time))
	fTimeProfile_mu[time]++;
      else
	fTimeProfile_mu[time] = 1;
      fNpe_mu += (1. / fSampFact);
    }
  }
}

void
HitStation::WriteADCBin (Int_t idpart, Int_t bin, Double_t bincontent)
{
  Short_t high;
  Short_t low;

  high = (Short_t) (bincontent * (DYNODETOANODE / VEMNORMALIZATION));
  low = (Short_t) (bincontent / VEMNORMALIZATION);

  fADC[1][bin] = low;
  fADC[0][bin] = high;
  fNpe_mu += (Int_t) fADC_mu[0][bin];

  if (theConfig ()->MuEmFlag) {
    if (idpart == 1 || idpart == -1) {
      fADC_em[1][bin] = low;
      fADC_em[0][bin] = high;
      fNpe_ph += fADC_em[0][bin];
    }
    if (idpart == 2 || idpart == -2) {
      fADC_em[1][bin] = low;
      fADC_em[0][bin] = high;
      fNpe_el += fADC_em[0][bin];
    }

    if (idpart == 3 || idpart == -3) {
      fADC_mu[1][bin] = low;
      fADC_mu[0][bin] = high;
      fNpe_mu += fADC_mu[0][bin];

    }
  }

}

void
HitStation::WriteDoubleADCBin (Int_t idpart, Int_t bin, Double_t bincontent)
{

  fDoubleADC[1][bin] += bincontent;
  fDoubleADC[0][bin] += bincontent * DYNODETOANODE;

  if (theConfig ()->MuEmFlag) {
    if (idpart == 1 || idpart == -1) {
      fDoubleADC_em[1][bin] += bincontent;
      fDoubleADC_em[0][bin] += bincontent * DYNODETOANODE;

    }
    if (idpart == 2 || idpart == -2) {
      fDoubleADC_em[1][bin] += bincontent;
      fDoubleADC_em[0][bin] += bincontent * DYNODETOANODE;

    }

    if (idpart == 3 || idpart == -3) {
      fDoubleADC_mu[1][bin] += bincontent;
      fDoubleADC_mu[0][bin] += bincontent * DYNODETOANODE;
    }
  }
}

void
HitStation::CountPart (Int_t idpart, Double_t wpart_det, Double_t time)
{


  fNPartTot += (Int_t) wpart_det;
  if (idpart == 1 || idpart == -1)
    fNph += (Int_t) wpart_det;

  else if (idpart == 2 || idpart == -2)
    fNel += (Int_t) wpart_det;

  else if (idpart == 3 || idpart == -3) {
    fNmu += (Int_t) wpart_det;
    fMuTimes.push_back ((Int_t) time);
  }

}

void
HitStation::WriteShortADCTraces ()
{
  Int_t npe = 0;

  for (Int_t i = 0; i < MAXNUMBEROFADCBINS; i++) {
    Short_t high, low, high_mu, low_mu, high_em, low_em;

    high = (Short_t) (fDoubleADC[0][i] / VEMNORMALIZATION);
    low = (Short_t) (fDoubleADC[1][i] / VEMNORMALIZATION);

    if (low > NBFADCCHANNEL)
      fADC[1][i] = NBFADCCHANNEL;
    else if (low < 0)
      fADC[1][i] = 0;
    else
      fADC[1][i] = low;

    if (high > NBFADCCHANNEL)
      fADC[0][i] = NBFADCCHANNEL;
    else if (high < 0)
      fADC[0][i] = 0;
    else
      fADC[0][i] = high;

    npe += fADC[0][i];

    if (theConfig ()->MuEmFlag) {
      high_mu = (Short_t) (fDoubleADC_mu[0][i] / VEMNORMALIZATION);
      low_mu = (Short_t) (fDoubleADC_mu[1][i] / VEMNORMALIZATION);
      high_em = (Short_t) (fDoubleADC_em[0][i] / VEMNORMALIZATION);
      low_em = (Short_t) (fDoubleADC_em[1][i] / VEMNORMALIZATION);

      if (low_mu > NBFADCCHANNEL)
	fADC_mu[1][i] = NBFADCCHANNEL;
      else if (low_mu < 0)
	fADC_mu[1][i] = 0;
      else
	fADC_mu[1][i] = low_mu;

      if (high_mu > NBFADCCHANNEL)
	fADC_mu[0][i] = NBFADCCHANNEL;
      else if (high_mu < 0)
	fADC_mu[0][i] = 0;
      else
	fADC_mu[0][i] = high_mu;

      if (low_em > NBFADCCHANNEL)
	fADC_em[1][i] = NBFADCCHANNEL;
      else if (low_em < 0)
	fADC_em[1][i] = 0;
      else
	fADC_em[1][i] = low_em;

      if (high_em > NBFADCCHANNEL)
	fADC_em[0][i] = NBFADCCHANNEL;
      else if (high_em < 0)
	fADC_em[0][i] = 0;
      else
	fADC_em[0][i] = high_em;
    }
  }

  fNpe = npe;

}

void
HitStation::ReadParticleInTank (ParticleInTank * p)
{
  Int_t nph = p->fNCherPhot;
  Int_t ndirect[NPM];

  for (int ipm = 0; ipm < NPM; ipm++)
    ndirect[ipm] = 0;
  for (Int_t np = 0; np < nph; np++) {
    Int_t ipm = p->fIpm[np];
    Int_t time = (Int_t) p->fTCherPhot[np];

    if (p->fNReflexions[np] == 0)
      ndirect[ipm - 1]++;
    //if(time<0 || ipm <=0) cout<<"hahaha"<< ipm<<" "<<time<<endl;
    if (ipm > 0) {
      fPMT[ipm - 1].WritePE (p->fId, time /*+TIMEOFFSET */ );
      WritePE (p->fId, time /*+TIMEOFFSET */ );

      if (time < fFirstPelBin) {
	fFirstPelBin = (int) time;
	//    firsttime[ipm-1] =(int) (time + time_plane);

	fFirstPelTime = (int) (time - Z_sf () / CLIGHT);

      }

      if (time < fPMT[ipm - 1].fFirstPelBin) {
	fPMT[ipm - 1].fFirstPelBin = (int) time;
	//    firsttime[ipm-1] =(int) (time + time_plane);

	fPMT[ipm - 1].fFirstPelTime = (int) (time - Z_sf () / CLIGHT);

      }
    }
  }

  int sumdirect = 0;
  for (int ipm = 0; ipm < NPM; ipm++) {
    fPMT[ipm].DirectLight (ndirect[ipm]);
    sumdirect += ndirect[ipm];
  }

  fNpe_direct += sumdirect;

}

void
HitStation::SumADCTraces ()
{

  for (Int_t ipm = 0; ipm < NPM; ipm++) {
    for (Int_t ig = 0; ig < 2; ig++) {
      for (Int_t ib = 0; ib < MAXNUMBEROFADCBINS; ib++) {
	fADC[ig][ib] += fPMT[ipm].fADC[ig][ib];
	if (theConfig ()->MuEmFlag) {
	  fADC_em[ig][ib] += fPMT[ipm].fADC_em[ig][ib];
	  fADC_mu[ig][ib] += fPMT[ipm].fADC_mu[ig][ib];
	}
	if (theConfig ()->ElecMode == "SHOWSAT" && ig == 1)
	  fADC_ns[ib] += fPMT[ipm].fADC_ns[ib];
      }
    }
  }

}

void
HitStation::DoDetSim (string simmode, double costheta_sh, double sintheta_sh,
		      double cosphi_sh, double sinphi_sh,
		      const SampStation & sst)
{
  Double_t eps = 0.001;
  Double_t sigt = 0.1;
  int npel_ph, npel_el, npel_mu;
  Int_t firsttime[NPM], firstbin[NPM];

  npel_ph = 0;
  npel_el = 0;
  npel_mu = 0;

  for (int ipm = 0; ipm < NPM; ipm++) {
    firsttime[ipm] = 0;
    firstbin[ipm] = MAXNUMBEROFTIMEBINS;
  }

  Int_t nparttot = sst.fNPartTot;
  if (nparttot == 0)
    return;

  Double_t sampfact;
  if (nparttot > NMAXPARTINTANK)
    sampfact = (Double_t) NMAXPARTINTANK / (Double_t) nparttot;
  else
    sampfact = 1.;
  this->WriteSamplingFactor (sampfact);

  Int_t npart = sst.fPartList.size ();


 Double_t maxdiftimepart=-99999999.99;
  Double_t mindiftimepart=99999999.99;

  for(Int_t np=0;np<npart;np++){
    Double_t timepart = sst.fPartList[np].fT;
    Double_t time_plane = sst.fPartList[np].fT_plane;
    if(timepart-time_plane<mindiftimepart)mindiftimepart=timepart-time_plane;
    if(timepart-time_plane>maxdiftimepart)maxdiftimepart=timepart-time_plane;
  }
  if(mindiftimepart>10000
     || (maxdiftimepart>10000 && maxdiftimepart-mindiftimepart<10000)){
    cerr<<"Correting Front time due to too large offset "
    <<mindiftimepart<<" "<<maxdiftimepart<<endl;
    mindiftimepart-=10;
  }
  else mindiftimepart=0;


  // loop on all particles in tank ( each particle has a weight )

  for (Int_t np = 0; np < npart; np++) {
    //      Particle &p= sst.fPartList[np];

    if (!np) {
      cerr.width (3);
      cerr << (int) (np * 1. / npart) << "%";
    } else {
      cerr << "\b\b\b\b";
      cerr.width (3);
      cerr << floor (np * 100. / ((npart - 1) * 1.)) << "%";
      fflush (stdout);
    }

    Int_t idpart = sst.fPartList[np].fId;
    Double_t cxpart = sst.fPartList[np].fUX;
    Double_t cypart = sst.fPartList[np].fUY;
    Double_t czpart = sst.fPartList[np].fUZ;
    Double_t epart = sst.fPartList[np].fE;
    Double_t timepart = sst.fPartList[np].fT;
    // Double_t time_plane = sst.fPartList[np].fT_plane;
   
    Double_t time_plane = sst.fPartList[np].fT_plane+mindiftimepart; //COBB
 // if(time_plane>timepart)cout<<"attention "<<timepart<<" "<<time_plane<<endl;
    Double_t timr0 = max (timepart - time_plane, 0.);
    Int_t n_entry_top = (Int_t) sst.fPartList[np].fWtop;
    Int_t n_entry_side = (Int_t) sst.fPartList[np].fWside;

    // the actual number of "clones" is bounded
    // if necessary, a weight is put
    Int_t nt = min (n_entry_top, NCLONES);
    Int_t ns = min (n_entry_side, NCLONES);

    //   top entries ---------------------------

    Double_t wpart_det;
    wpart_det = (Double_t) n_entry_top *sampfact / (Double_t) nt;

    for (Int_t n = 0; n < nt; n++) {
      //  random entry point
      Double_t ri = (STATION_RADIUS - eps) * sqrt (Rand ());
      Double_t azi = 2 * PI * Rand ();
      Double_t xx = ri * cos (azi);
      Double_t yy = ri * sin (azi);
      Double_t zz = STATION_HEIGHT - eps;
      //   smoothed time (if multiple entries from one input particle)
      Double_t timr = timr0;
      if (nt + ns > 1)
	timr = timr0 * exp (sigt * Rand ("GAUSS", 0, 1));
      //  delay of entry related to the position on the wall (time=0 at origin)
      timr = timr - (sintheta_sh * (xx * cosphi_sh + yy * sinphi_sh)
		     - costheta_sh * zz) / CLIGHT;

      ParticleInTank *pclone = new ParticleInTank (fId, idpart, xx, yy, zz,
						   cxpart, cypart, czpart,
						   timr, epart, wpart_det);

      this->CountPart (idpart, wpart_det / sampfact, timr);
      if (simmode != "SAMPLE")
	pclone->DoDetSim (simmode);
      //if(simmode != "SAMPLE") pclone->DoDetSim();

      if (simmode == "DETAILED") {
	ReadParticleInTank (pclone);

	//here we count the pel for a summary
	Int_t nph = pclone->fNCherPhot;

	if (abs (idpart) == 1)
	  npel_ph += nph;
	else if (abs (idpart) == 2)
	  npel_el += nph;
	else if (abs (idpart) == 3)
	  npel_mu += nph;

      }

      else			//fast simulation
      {
	Int_t nbins = pclone->fBinContent.size ();
	for (Int_t nb = 0; nb < nbins; nb++) {
	  Double_t bincontent = pclone->fBinContent[nb];
	  this->WriteDoubleADCBin (idpart, nb, bincontent);
	}
      }
      delete pclone;

    }				//end of loop on top entries

    //   side entries ---------------------------

    wpart_det = (Double_t) n_entry_side *sampfact / (Double_t) ns;
    Double_t cxy = sqrt (cxpart * cxpart + cypart * cypart);


    for (Int_t n = 0; n < ns; n++) {
      //  random entry point
      Double_t yy = (STATION_RADIUS - eps) * (2 * Rand () - 1.);
      Double_t xxx =
	sqrt ((STATION_RADIUS - eps) * (STATION_RADIUS - eps) - yy * yy);
      Double_t xx = (-xxx * cxpart + yy * cypart) / cxy;
      yy = (-xxx * cypart - yy * cxpart) / cxy;
      Double_t zz = (STATION_HEIGHT - eps) * Rand ();

      //   smoothed time (if multiple entries from one input particle)
      Double_t timr = timr0;
      if (nt + ns > 1)
	timr = timr0 * exp (sigt * Rand ("GAUSS", 0, 1));

      //  delay of entry related to the position on the wall (time=0 at origin)
      timr = timr - (sintheta_sh * (xx * cosphi_sh + yy * sinphi_sh)
		     - costheta_sh * zz) / CLIGHT;

      ParticleInTank *pclone =
	new ParticleInTank (fId, idpart, xx, yy, zz, cxpart, cypart, czpart,
			    timr, epart, wpart_det);

      this->CountPart (idpart, wpart_det / sampfact, timr);

      if (simmode != "SAMPLE")
	pclone->DoDetSim (simmode);
      //if(simmode != "SAMPLE") pclone->DoDetSim();

      if (simmode == "DETAILED") {
	ReadParticleInTank (pclone);
	//here we count the pel for a summary
	Int_t nph = pclone->fNCherPhot;

	if (abs (idpart) == 1)
	  npel_ph += nph;
	else if (abs (idpart) == 2)
	  npel_el += nph;
	else if (abs (idpart) == 3)
	  npel_mu += nph;
      } else			//fast simulation
      {
	Int_t nbins = pclone->fBinContent.size ();
	for (Int_t nb = 0; nb < nbins; nb++) {
	  Double_t bincontent = pclone->fBinContent[nb];
	  this->WriteDoubleADCBin (idpart, nb, bincontent);
	}
      }

      delete pclone;
    }				//end of loop on side entries

  }				//end of loop on part in a tank


}

void
HitStation::SetTime0 ()
{
  int firstADCbin = MAXNUMBEROFTIMEBINS;
  int firstpm = 0;

  for (int ipm = 0; ipm < NPM; ipm++) {
    if (firstADCbin > fPMT[ipm].fFirstADCBin) {
      firstADCbin = fPMT[ipm].fFirstADCBin;
      firstpm = ipm;
    }
  }
  //fTime0=fPMT[firstpm].fFirstPelTime;
  // fTime0=fPMT[firstpm].fFirstPelTime;//  + SAMPLINGTIME * firstADCbin -fPMT[firstpm].fFirstPelBin;
  fTime0 = fFirstPelTime + SAMPLINGTIME * fTriggerBin - fFirstPelBin;
}

void
HitStation::SetTime0 (double time0)
{
  fTime0 = time0;
}

void
HitStation::DoElecSim ()
{

  for (Int_t ipm = 0; ipm < NPM; ipm++) {
    fPMT[ipm].DoElecSim ();

  }
  SumADCTraces ();
}

int
HitStation::DoTrigSim ()
{
  int bintrig = 0;
  int t2trig = 0;
  if (TestToT ()) {
    bintrig = TestToT ();
    fT2ToT = 1;
    t2trig = 1;
    cerr << "bin ToT " << bintrig << endl;
  }

  if (TestThreshold (2)) {
    bintrig = TestThreshold (2);
    fT2Threshold = 1;
    fT1Threshold = 1;
    t2trig = 1;
    cerr << "bin T2Th " << bintrig << endl;
  }

  if (!t2trig)
    if (TestThreshold (1)) {
      bintrig = TestThreshold (1);
      fT1Threshold = 1;
      cerr << "bin T1Th " << bintrig << endl;

    }


  fTriggerBin = bintrig;

  return t2trig;
}

int
HitStation::TestToT ()
{

  int nmul = 0;
  int nbintrig;
  int nbintrigmax = 0;
  int flagtrig = 0;
  int trigbin = 0;

  for (int i = 0; i < MAXNUMBEROFADCBINS - TOTTRIGWINDOWSIZE; i++) {
    nbintrig = 0;
    for (int j = 0; j < TOTTRIGWINDOWSIZE; j++) {
      nmul = 0;
      for (Int_t ipm = 0; ipm < NPM; ipm++) {
	double vem = fPMT[ipm].fADC[0][i + j] / VEMPEAKVALUEINADC;
	if (vem > TOTTRIGTHRESHOLD)
	  nmul++;
      }
      if (nmul >= TOTMULTIPLICITY)
	nbintrig++;

      //  cout<<i+j<<" "<<nmul<<" "<<nbintrig<<endl;

      if (nbintrig > nbintrigmax) {
	nbintrigmax = nbintrig;
	if (nbintrigmax > TOTTRIGNUMBEROFBINS && !flagtrig) {
	  flagtrig = 1;
	  trigbin = i + j;
	  //  cout<<"attention"<<trigbin<<endl;
	}
      }

    }
  }				// loop on time slice of the Tot Window

  fNbOfBinsOverThres = nbintrigmax;
  return (trigbin);


}

int
HitStation::TestThreshold (int trigtype)
{
  int flagtrig = 0;
  double threshold = 0;
  int bintrig = 0;
  int nmul = 0;
  if (trigtype == 1)
    threshold = T1TRIGTHRESHOLD;
  if (trigtype == 2)
    threshold = T2TRIGTHRESHOLD;
  for (int i = 0; i < MAXNUMBEROFADCBINS; i++) {
    nmul = 0;
    for (Int_t ipm = 0; ipm < NPM; ipm++) {
      double vem = fPMT[ipm].fADC[0][i] / VEMPEAKVALUEINADC;

      if (vem > threshold)
	nmul++;
      // cout<<"mul "<<i<<" "<<vem<<" "<<nmul<<endl;

    }
    if (nmul >= THRESHOLDMULTIPLICITY) {
      flagtrig = 1;
      bintrig = i;
      break;
    }
  }

  return (bintrig);

}


double
HitStation::SignalInVem ()
{
  double signal = 0.;
  int satur = 0;
  for (int i = 0; i < MAXNUMBEROFADCBINS; i++) {
    signal += (int) fADC[0][i];
    if (fADC[0][i] >= 1023) {
      satur = 1;
      break;
    }
  }

  if (satur == 1) {
    signal = 0;
    for (int i = 0; i < MAXNUMBEROFADCBINS; i++)
      signal += (int) fADC[1][i];
    signal *= DYNODETOANODE;
  }

  fSignalInVem = signal / (VEMCHARGEVALUEINADC * (float)NPM);

  return fSignalInVem;
}
