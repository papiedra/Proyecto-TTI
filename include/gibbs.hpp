#ifndef _GIBBS_
#define _GIBBS_
/**
 * @file gibbs.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera gibbs
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
 * @return tuple<Matrix&, double, double, double,string> 
 */
tuple<Matrix&, double, double, double,string> gibbs(Matrix r1, Matrix r2, Matrix r3);

#endif