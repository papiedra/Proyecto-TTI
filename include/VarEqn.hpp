#ifndef _VAREQN_
#define _VAREQN_

#include "..\include\matrix.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMat.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"


Matrix& VarEqn(double x, Matrix yPhi);

#endif