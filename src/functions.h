 /***************************************************
  * H. G. Asorey                                    *
  * asoreyh at cab.cnea.gov.ar                      *
  ***************************************************/

#include "ran2.cpp"

//Declaraciones
#ifdef CUBE
inline bool isInCube(vect mcar);
#endif
inline bool isIn(vect mcyl);
double pdf(double x);
double AcceptRej(double (*function)(double x), double max);
int getType(vect mcyl, bool in);
int getTypeCube(vect mcar, bool in);
vect cyl2car (vect cyl);
vect car2cyl (vect car);
vect cyl2sph (vect cyl);
vect sph2cyl (vect sph);
vect sph2car (vect sph);
vect car2sph (vect car);

//Funciones
vect car2sph (vect car) {
  vect sph;
  sph.r1 = sqrt(car.r1 * car.r1 + car.r2 * car.r2 + car.r3 * car.r3);
  sph.r2 = (sph.r1 > 0. ? acos(car.r3 / sph.r1) : M_PI / 2.);
  sph.r3 = (car.r1 > 0. ? atan(car.r2 / car.r1) : M_PI / 2.);
  return sph;
}

vect sph2car (vect sph) {
  vect car;
  car.r1 = sph.r1 * sin(sph.r2) * cos(sph.r3);
  car.r2 = sph.r1 * sin(sph.r2) * sin(sph.r3);
  car.r3 = sph.r1 * cos(sph.r2);
  return car;
}

vect sph2cyl (vect sph) {
  vect cyl;
  cyl.r1 = sph.r1 * sin(sph.r2);
  cyl.r2 = sph.r3;
  cyl.r3 = sph.r1 * cos(sph.r2);
  return cyl;
}

vect cyl2sph (vect cyl) {
  vect sph;
  sph.r1 = sqrt(cyl.r1 * cyl.r1 + cyl.r3 * cyl.r3);
  sph.r2 = (cyl.r3 > 0. ? atan(cyl.r1 / cyl.r3) : M_PI / 2.);
  sph.r3 = cyl.r2;
  return cyl;
}


vect car2cyl (vect car) {
  vect cyl;
  cyl.r1 = sqrt(car.r1 * car.r1 + car.r2 * car.r2);
  if (car.r1 == 0. && car.r2 == 0.) cyl.r2 = 0.;
  else if (car.r1 >= 0) cyl.r2 = asin(car.r2 / cyl.r1);
  else if (car.r1 <0) cyl.r2 = M_PI - asin(car.r2 / cyl.r1);
  cyl.r3 = car.r3;
  return cyl;
}

vect cyl2car (vect cyl) {
  vect car;
  car.r1 = cyl.r1 * cos(cyl.r2);
  car.r2 = cyl.r1 * sin(cyl.r2);
  car.r3 = cyl.r3;
  return car;
}

double pdf(double x) {
//  return (cos(x)*cos(x));
  return (cos(x)*cos(x)*sin(x));
//  return (pow(cos(x),3)*sin(x));
//  return (cos(x)*cos(x)*cos(x)*cos(x)*sin(x));
}

double AcceptRej(double (*function)(double x), double max){
  double y,u;
  do {
    y = ran2(&seedvy) * max;
    u = ran2(&seedvu);
  } while (u > (*function)(y));
  return y; 
}

inline bool isIn(vect mcyl) {
  if (mcyl.r3 > 0. && mcyl.r3 < tankH)
    if (mcyl.r1 < tankR)
      return true;
  return false;
}

#ifdef CUBE
inline bool isInCube(vect mcar) {
  if (mcar.r3 > 0. && mcar.r3 <= tankH)
    if (mcar.r2 >= -tankB2 && mcar.r2 <= tankB2)
      if (mcar.r1 >= -tankA2 && mcar.r1 <= tankA2)
        return true;
  return false;
}
#endif

int getType(vect mcyl, bool in) {
  if (!in) { // the muon is entering the tank
    if (mcyl.r1 < tankR) return 1; // the muon is entering by the top
    else return 2; // the muon is entering by side
  }
  else { // the muon is going out
    if ( type == 1 ) {
      if (mcyl.r1 < tankR) return 1; // top-bottom
      else return 2; //top-side
    }
    if ( type == 2 && mcyl.r1 < tankR) return 3; // side-bottom
    else return 4; // side-side
  }
}

#ifdef CUBE
int getTypeCube(vect mcar, bool in) {
  if (!in) { // the muon is entering the tank
//  if (mcar.r1 > -tankA2 && mcar.r1 < tankA2 && mcar.r2 > -tankB2 && mcar.r2 < tankB2 )
    if (mcar.r3 >= (tankH-1))
      return 1; // the muon is entering by the top
    else 
      return 2; // the muon is entering by side
  }
  else { // the muon is going out
    if ( type == 1 ) {
      if (mcar.r1 > -tankA2 && mcar.r1 < tankA2 && mcar.r2 > -tankB2 && mcar.r2 < tankB2 ) 
        return 1; // top-bottom
      else 
        return 2; //top-side
    }
    if ( type == 2 && mcar.r1 > -tankA2 && mcar.r1 < tankA2 && mcar.r2 > -tankB2 && mcar.r2 < tankB2 )
      return 3; // side-bottom
    else 
      return 4; // side-side
  }
}
#endif
