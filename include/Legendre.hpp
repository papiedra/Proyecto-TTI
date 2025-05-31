#ifndef _LEGENDRE_
#define _LEGENDRE_
/**
 * @file Legendre.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera legendre
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "matrix.hpp"
#include <tuple>
/**
 * @brief 
 * 
 * @param n 
 * @param m 
 * @param fi 
 * @return std::tuple<Matrix, Matrix> 
 */
std::tuple<Matrix, Matrix> legendre(int n, int m, double fi);

#endif