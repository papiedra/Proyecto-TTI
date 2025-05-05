#ifndef _MEASUPDATE_
#define _MEASUPDATE_
#include <cmath>
#include <tuple>
#include "..\include\matrix.hpp"

std::tuple<Matrix, Matrix, Matrix> MeasUpdate(Matrix& x,Matrix& z,Matrix& g,Matrix& s,Matrix& G,Matrix& P,double n);

#endif