#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_
/**
 * @file JPL_Eph_DE430.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera JPL_Eph_DE430
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <cmath>
#include <tuple>
#include "..\include\matrix.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\global.hpp"

/**
 * @brief 
 * 
 * @param Mjd_TDB 
 * @return std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> 
 */
std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);


#endif