#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_
/**
 * @file AccelHarmonic.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera AccelHarmonic
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "matrix.hpp"
#include <cmath>
/**
 * @brief 
 * 
 * @param r 
 * @param E 
 * @param n_max 
 * @param m_max 
 * @return Matrix& 
 */
Matrix& AccelHarmonic(Matrix& r, Matrix E, int n_max,int m_max);



#endif