#ifndef _HGIBBS_
#define _HGIBBS_
/**
 * @file hgibbs.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera hgibbs
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\unit.hpp"
#include <cmath>
#include <tuple>
/**
 * @brief 
 * 
 * @param r1 
 * @param r2 
 * @param r3 
 * @param Mjd1 
 * @param Mjd2 
 * @param Mjd3 
 * @return tuple<Matrix&, double, double, double, string> 
 */
tuple<Matrix&, double, double, double, string> hgibbs(Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3);

#endif