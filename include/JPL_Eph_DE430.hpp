#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include <cmath>
#include <tuple>
#include "..\include\matrix.hpp"

std::tuple<Matrix ,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix> JPL_Eph_DE430(double Mjd_TDB);


#endif