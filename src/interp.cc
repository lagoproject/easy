/* 3D interpolation code */
/* la idea es trabajar con álgebra vectorial, extendiendo lo que se hace en 2D a un vector de dimension n */


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
const double units    =  10.;  // in milimeters

// vamos a encontrar por interpolación el punto donde una línea intersecta a un cilindro de radio

const double tankR    = 1784. / units;
const double tankH    = 1200. / units;
const double dt       = 1.; // temporal step 

struct vect {
  double r1;
  double r2;
  double r3;
} coordinate;

#include "interp_func.h"

int main() {
/*
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
*/


  /* handlers */
  vect pcar = { 0., 0., 0.};  // photon position in cartesian
  vect pcyl = { 0., 0., 0.};  // photon position in cylindrical

  vect ppcar = { 0., 0., 0.};  // photon position in cylindrical
  vect ppcyl = { 0., 0., 0.};  // photon position in cylindrical

  vect vcar = { 0., 0., 0.};  // photon speed in cartesian
  vect vsph = { 0., 0., 0.};  // photon speed in spherical

  //photon at the top
  pcar.r1 = 0.;
  pcar.r2 = 0.;
  pcar.r3= tankH;
  pcyl = car2cyl(pcar);

  // let's use 70 going down at phi=0, |v| = 1.; in sphericals: 
  vsph.r1 = -1.;
  vsph.r2= 70.*deg2rad;
  vsph.r3= 20.*deg2rad;
  vcar = sph2car(vsph);
 
  do {
    // move the photon 
    ppcyl = pcyl;
    pcar.r1 += vcar.r1 * dt;
    pcar.r2 += vcar.r2 * dt;
    pcar.r3 += vcar.r3 * dt;
    pcyl = car2cyl(pcar);
  } while (isIn(pcyl));
  ppcar = cyl2car(pcyl);
  vect ptnk=intersect(ppcar,pcar); // intersección con el tanque
  ptnk = car2cyl(ptnk);
  cout 
    << ppcyl.r1 << " "
    << ppcyl.r2 << " "
    << ppcyl.r3 << " " << endl 
    << ptnk.r1 << " "
    << ptnk.r2 << " "
    << ptnk.r3 << " " << endl 
    << pcyl.r1 << " "
    << pcyl.r2 << " "
    << pcyl.r3 << " " << endl;

  // pcyl tiene la ultima posición (fuera), ppcyl tiene la ultima posición dentro


} //main
