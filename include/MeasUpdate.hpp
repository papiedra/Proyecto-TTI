#ifndef _MEASUPDATE_
#define _MEASUPDATE_
/**
 * @file MeasUpdate.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera MeasUpdate
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include <tuple>
/**
 * @brief 
 * 
 * @param x 
 * @param z 
 * @param g 
 * @param s 
 * @param G 
 * @param P 
 * @param n 
 * @return tuple<Matrix&, Matrix&, Matrix&> 
 */
tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n);

#endif