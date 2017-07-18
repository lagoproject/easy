
/*
  Utils.cc

  implementation file for Utils.h

*/


#include "Utils.h"
#include "TROOT.h"
#include "TRandom3.h"


void UTMToXY(double north,double east,double *x,double *y)
{
  static const double RefNorthing = 6.08276E+06;
  static const double RefEasting = 459630.;
  static const double Alpha = 2.52E-03;
  static const double Beta = 6.03E-04;

  *x = ((1+Beta)*(east-RefEasting)+Alpha*(north-RefNorthing));
  *y = ((1+Beta)*(north-RefNorthing)-Alpha*(east-RefEasting));
}

double Rand(string mode, double mean, double sigma) 
{
  double ret=0.;
  static TRandom3 *_rand = NULL;
  if (!_rand)
    {
      _rand = new TRandom3();
      _rand->SetSeed(0);
    }

  if (mode == "")
    {
      ret = _rand->Rndm();
    }
  else if (mode == "POISSON")
    {
      ret = _rand->Poisson(mean);
    }
  else if (mode == "GAUSS")
    {
      ret = _rand->Gaus(mean, sigma);
    }

  return ret;
}


