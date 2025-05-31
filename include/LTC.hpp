#ifndef _LTC_
#define _LTC_
/**
 * @file LTC.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera LTC
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
 * @param lon 
 * @param lat 
 * @return Matrix 
 */
Matrix LTC(double lon, double lat);

#endif