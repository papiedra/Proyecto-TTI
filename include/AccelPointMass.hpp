#ifndef _ACCELPOINTMASS_
#define _ACCELPOINTMASS_
/**
 * @file AccelPointMass.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera AccelPointMass
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
 * @param r 
 * @param s 
 * @param GM 
 * @return Matrix& 
 */
Matrix& AccelPointMass(Matrix &r, Matrix &s, double GM);

#endif