#ifndef _VAREQN_
#define _VAREQN_
/**
 * @file VarEqn.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera VarEqn
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
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

/**
 * @brief 
 * 
 * @param x 
 * @param yPhi 
 * @return Matrix& 
 */
Matrix& VarEqn(double x, Matrix yPhi);

#endif