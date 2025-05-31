#ifndef _NUTANGLES_
#define _NUTANGLES_
/**
 * @file NutAngles.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera NutAngles
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <tuple>
#include "SAT_Const.hpp"
#include "matrix.hpp"
#include <cmath>

/**
 * @brief 
 * 
 * @param Mjd_TT 
 * @return std::tuple<double, double> 
 */
std::tuple<double, double> NutAngles(double Mjd_TT);

#endif