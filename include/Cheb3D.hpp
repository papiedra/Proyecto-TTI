#ifndef _Cheb3D_
#define _Cheb3D_
/**
 * @file Cheb3D.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera Cheb3D
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include <cmath>
/**
 * @brief 
 * 
 * @param t 
 * @param N 
 * @param Ta 
 * @param Tb 
 * @param Cx 
 * @param Cy 
 * @param Cz 
 * @return Matrix& 
 */
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz);

#endif