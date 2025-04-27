#ifndef _IERS_
#define _IERS_
#include <tuple>
#include "..\include\matrix.hpp"

std::tuple<double,  double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp);

#endif