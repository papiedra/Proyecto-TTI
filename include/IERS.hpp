#ifndef _IERS_
#define _IERS_
/**
 * @file IERS.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera IERS
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <tuple>
#include "..\include\matrix.hpp"
/**
 * @brief 
 * 
 * @param eop 
 * @param Mjd_UTC 
 * @param interp 
 * @return std::tuple<double,  double, double, double, double, double, double, double, double> 
 */
std::tuple<double,  double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp);

#endif