#ifndef _GEODETIC_
#define _GEODETIC_
/**
 * @file Geodetic.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera Geodetic
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <tuple>
/**
 * @brief 
 * 
 * @param r 
 * @return tuple<double, double, double> 
 */
tuple<double, double, double> Geodetic(Matrix r);

#endif