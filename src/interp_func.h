 /***************************************************
  * H. G. Asorey                                    *
  * asoreyh at cab.cnea.gov.ar                      *
  ***************************************************/

//Declaraciones
inline bool isIn(vect mcyl);
vect cyl2car (vect cyl);
vect car2cyl (vect car);
vect cyl2sph (vect cyl);
vect sph2cyl (vect sph);
vect sph2car (vect sph);
vect car2sph (vect car);
int out(vect rout);//rout en cilindricas
vect intersect (vect r1, vect r2);
inline double mod(vect r);
double dist(vect p1,vect p2);

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
  return sph;
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

inline double mod(vect r) {
  return (sqrt(r.r1*r.r1+r.r2*r.r2+r.r3*r.r3));
}

double dist(vect p1,vect p2) {
  vect d;
  d.r1=p2.r1-p1.r1;
  d.r2=p2.r2-p1.r2;
  d.r3=p2.r3-p1.r3;
  return (mod(d));
}

inline bool isIn(vect mcyl) {
  if (mcyl.r3 > 0. && mcyl.r3 < tankH)
    if (mcyl.r1 < tankR)
      return true;
  return false;
}

int out(vect rout) { //rout en cilindricas
/* determina si salio por el fondo (-1), techo (1) o lateral (0) */
  if (rout.r1 <= tankR) {
    if (rout.r3 > tankH) // techo
      return 1;
    else // fondo
      return -1;
  }
  return 0; //dado que esta afuera, y no salio por arriba ni por abajo, salio por el lateral
}

vect intersect (vect p1, vect p2) { //p1=ultimo dentro y p2=primero afuera en cartesianas
  double dx=p1.r1-p2.r1;
  double dy=p1.r2-p2.r2;
  double dz=p2.r3-p1.r3;
  double xy, sx, sy, sd;
  double dr12, dr1a, dr1b;
  vect ptnka, ptnkb, ptnk;
  switch (out(car2cyl(p2))) {
    case 0: // lateral, intersección de un cilindro con una esfera, dos posibles puntos, a y b.
      xy = p2.r1*p1.r2-p1.r1*p2.r2;
      sd = dx*dx+dy*dy;
      sx = sqrt(dx*dx*(tankR*tankR*sd-xy*xy));
      sy = p1.r1*p1.r1*p1.r1*p2.r2-p2.r1*p2.r1*p2.r1*p1.r2+p1.r1*p2.r1*p2.r1*(2.*p1.r2+p2.r2)-p1.r1*p1.r1*p2.r1*(p1.r2+2.*p2.r2);
      ptnka.r1 = (dy*xy+sx)/sd;
      ptnkb.r1 = (dy*xy-sx)/sd;
      ptnka.r2 = (sy+sx*dy)/(dx*sd);
      ptnkb.r2 = (sy-sx*dy)/(dx*sd);
      dr12 = sqrt(sd);
      dr1a = sqrt((p1.r1-ptnka.r1)*(p1.r1-ptnka.r1)+(p1.r2-ptnka.r2)*(p1.r2-ptnka.r2));
      dr1b = sqrt((p1.r1-ptnkb.r1)*(p1.r1-ptnkb.r1)+(p1.r2-ptnkb.r2)*(p1.r2-ptnkb.r2));
      ptnka.r3 = p1.r3 + dz * dr1a / dr12;
      ptnkb.r3 = p1.r3 + dz * dr1b / dr12;
      if (dist(ptnka,p2) <= dist(ptnkb,p2))
        ptnk=ptnka;
      else
        ptnk=ptnkb;
      break;
    case 1: // techo
      break;
    case 2:
      break;
  }
  return ptnk;
}
