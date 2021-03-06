 /***************************************************
  * H. G. Asorey - 2011                             *
  * asoreyh at gmail.com                            *
  *                                                 *
  * determines the trace lenght for particles with  *
  * cos^2(q) distribution                           *
  ***************************************************/

/****************************************************
 * History log                                      *
 *                                                  *
 * v0r0 - 09 Oct 2011:                              *
 * first version - tag 152 on svn                   *
 *                                                  *
 * v0r1 - 17 Oct 2011:                              *
 * previous version does not represent correctly    *
 * the traces since the flux and not the speed is   *
 * modulated by cos^2(q). In this version the speed *
 * is uniform and the trace is weighted by cos^2(q) *
 *                                                  *
 * v1r0 - Ready. Next step: read external file for  *
 * particle injection                               *
 *                                                  *
 *                                                  *
 *                                                  *
 *                                                  *
 ****************************************************/

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

// config
const double deg2rad  = M_PI / 180.;
const double rad2deg  = 180. / M_PI;
const bool   internal = true; // use internal generator  or read muons from external file
const bool   extSeeds = true;  // use rand() to generate the seeds or just set them as -1
const double muons    =  1e6;  // numbers of showers to generate if internal is enabled
const double units    =  10.;  // in milimeters

// Parameters
/*
// nahuelito
const double tankR =  740. / units;
const double tankH = 1560. / units;
#define MU "nahuel.muo" 
#define TR "nahuel.trc"
#define DI "nahuel.dst"
*/

// boyita
/*
const double tankR =  800. / units;
const double tankH =  810. / units;
#define MU "boya.muo" 
#define TR "boya.trc"
#define DI "boya.dst"
*/

// auger
/*
const double tankR = 1784. / units;
const double tankH = 1200. / units;
#define MU "auger1.muo" 
#define TR "auger1.trc"
#define DI "auger1.dst"
*/

// guane-3
const double tankR = 515 / units; // detector radius
const double tankH = 500. / units; // detector height
#define MU "guane.muo" // muon crossing the detector distribution file
#define TR "guane.trc" // traces distribution file
#define DI "guane.dst" // all muons speed distribution file

// antartico
/*
const double tankR =  485. / units;
const double tankH = 1000. / units;
#define MU "antar.muo" 
#define TR "antar.trc"
#define DI "antar.dst"
*/

const double maxT  =   90. * deg2rad;
const int    maxL  = int(sqrt(4. * tankR * tankR + tankH * tankH)) + 1;
const int    types =    5; // type: 0: all; 1: top-bot; 2: top-side; 3:sid-bot; 4: side-side
const double skyR  =  10.0 * tankR;
const double skyH  =   2.2 * tankH;

// Random seeds

long int
  seedmx = -1,
  seedmy = -1,
  seedmz = -1,
  seedvy = -1,
  seedvu = -1,
  seedvf = -1;

struct vect {
  double r1;
  double r2;
  double r3;
} coordinate;

int type = 0;

#include "functions.h"

int main() {
  srand(time(NULL));
  if (extSeeds) {
    seedmx *= rand();
    seedmy *= rand();
    seedmz *= rand();
    seedvy *= rand();
    seedvu *= rand();
    seedvf *= rand();
  }

  ofstream dis,out,muo;
  char nmuo[256], nout[256], ndis[256]; 
  snprintf(nmuo,256, MU);
  snprintf(nout,256, TR);
  snprintf(ndis,256, DI);
  muo.open(nmuo);
  dis.open(ndis);

  vector < vector <double> > traces;
  traces.resize(maxL);
  for (int i = 0; i < maxL; i++)
    traces[i].resize(types);
  for (int i = 0; i < maxL; i++)
    for (int j = 0; j < types; j++)
      traces[i][j] = 0.;

  double
    c2    = 0.0, 
    trace = 0.0;

  int 
    cnt = 0,
    tot = 0;

  bool 
    in = false;

  vect mcar = { 0., 0., 0.};  // muon position in cartesian
  vect mcyl = { 0., 0., 0.};  // muon position in cylindrical
  vect pcyl = { 0., 0., 0.};  // muon position in cylindrical
  vect msph = { 0., 0., 0.};  // muon position in spherical
  vect ment = { 0., 0., 0.};  // muon position entering the tank

  vect vcar = { 0., 0., 0.};  // muon speed in cartesian
  vect vcyl = { 0., 0., 0.};  // muon speed in cylindrical
  vect vsph = { 0., 0., 0.};  // muon speed in spherical

  while (cnt < muons) {
    do {
      tot++;
      // put a muon in the sky
      mcar.r1 = (2. * ran2(&seedmx) - 1.) * skyR;
      mcar.r2 = (2. * ran2(&seedmy) - 1.) * skyR;
      mcar.r3 = skyH;
      mcyl = car2cyl(mcar);

      // and now determine the speed, 
      vsph.r1 = -1.;
      vsph.r2 = AcceptRej(pdf, maxT);
      vsph.r3 = ran2(&seedvf) * 2. * M_PI;
      vcar = sph2car(vsph);
      dis
        << tot     << " "
        << vsph.r1 << " "
        << vsph.r2 << " "
        << vsph.r3 << " "
        << endl;
 
      in = false;
      type = 1;
      trace = 0.;
      while (mcar.r3 >= 0.) {
        if (isIn(mcyl)) { 
          trace++;
          if (!in) {
            type = getType(pcyl, in);
            ment = cyl2car(pcyl);
            in = true;
          }
        }
        else if (in) {
          break;
        }
        // move the muon 
        pcyl = mcyl;
        mcar.r1 += vcar.r1;
        mcar.r2 += vcar.r2;
        mcar.r3 += vcar.r3;
        mcyl = car2cyl(mcar);
      }
    } while (trace <= 0.);
    trace = int(trace);
    type = getType(mcyl, in);
    muo
      << tot     << " "
      << ment.r1 << " "
      << ment.r2 << " "
      << ment.r3 << " "
      << mcar.r1 << " "
      << mcar.r2 << " "
      << mcar.r3 << " "
      << type    << " "
      << trace   << " "
      << endl;

    traces[trace][0]++;
    traces[trace][type]++;
    cnt++;
    if (!(cnt%(int(muons/100.)))) {
      cerr << cnt << "/" << tot << "\t" << int(100.*cnt/muons) << " %\r";
      out.open(nout);
      for (int i = 0; i < maxL; i++) 
        out 
          << i << " " // trace [cm] 
          << traces[i][0] << " " //total
          << traces[i][1] << " "
          << traces[i][2] << " "
          << traces[i][3] << " "
          << traces[i][4] << " "
          << endl;
      out.close();
    }
  }
  muo.close();
  dis.close();
  cerr << " particles : " << cnt << "/" << tot << endl;
} //main
