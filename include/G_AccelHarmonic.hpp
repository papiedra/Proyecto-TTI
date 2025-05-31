#ifndef _G_ACCELHARMONIC_
#define _G_ACCELHARMONIC_
/**
 * @file G_AccelHarmonic.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera AccelHarmonic
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <cmath>

#include "..\include\matrix.hpp"
#include "..\include\AccelHarmonic.hpp"
/**
 * @brief 
 * 
 * @param r 
 * @param U 
 * @param n_max 
 * @param m_max 
 * @return Matrix 
 */
Matrix G_AccelHarmonic(Matrix &r,Matrix &U,double n_max,double m_max);

#endif