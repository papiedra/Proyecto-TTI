#ifndef _ACCEL_
#define _ACCEL_

#include "..\include\matrix.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMat.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Mjday_TBD.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"


Matrix& Accel(double x, Matrix Y);

#endif