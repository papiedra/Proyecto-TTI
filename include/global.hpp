#ifndef _GLOBAL_
#define _GLOBAL_

#include "matrix.hpp"
#include <cmath>

typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern Param AuxParam;

void eop19620101(int c);
void GGM03S();
void DE430Coeff();
void auxparam();

#endif