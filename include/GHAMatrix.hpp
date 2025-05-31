#ifndef _GHAMATRIX_
#define _GHAMATRIX_
/**
 * @file GHAMatrix.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera GHAMatrix
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <cmath>

#include "..\include\matrix.hpp"
#include "..\include\R_z.hpp"
#include "..\include\gast.hpp"
/**
 * @brief 
 * 
 * @param Mjd_UT1 
 * @return Matrix 
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif