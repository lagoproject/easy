//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File Utils.h 

 */
//----------------------------------------------------------------------

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
using namespace std;

#include <string>

void UTMToXY(double north,double east,double *x,double *y);
double Rand(string mode = "", double mean = 0, double sigma = 1);

#endif
