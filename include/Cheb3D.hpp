#ifndef _Cheb3D_
#define _Cheb3D_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz);

#endif