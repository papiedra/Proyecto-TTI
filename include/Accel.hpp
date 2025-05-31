#ifndef _ACCEL_
#define _ACCEL_
/**
 * @file Accel.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera Accel
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
#include "..\include\Mjday_TBD.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
/**
 * @brief 
 * 
 * @param x 
 * @param Y 
 * @return Matrix& 
 */

Matrix& Accel(double x, Matrix Y);

#endif