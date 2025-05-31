#ifndef _AZELPA_
#define _AZELPA_
/**
 * @file AzElPa.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera AzElPa
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
 * @param s 
 * @return std::tuple<double,  double, Matrix, Matrix> 
 */
std::tuple<double,  double, Matrix, Matrix> AzElPa(Matrix& s);

#endif