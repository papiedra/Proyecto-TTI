#ifndef _DEINTEG_
#define _DEINTEG_
/**
 * @file DEInteg.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera DEInteg
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\sign_.hpp"

#include <cmath>
/**
 * @brief 
 * 
 * @param f 
 * @param t 
 * @param tout 
 * @param relerr 
 * @param abserr 
 * @param n_eqn 
 * @param y 
 * @return Matrix& 
 */
Matrix& DEInteg(Matrix& f(double t, Matrix y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);

#endif