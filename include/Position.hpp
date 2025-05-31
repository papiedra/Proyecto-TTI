#ifndef _POSITION_
#define _POSITION_
/**
 * @file Position.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera Position
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
/**
 * @brief 
 * 
 * @param lon 
 * @param lat 
 * @param h 
 * @return Matrix 
 */
Matrix Position(double lon, double lat, double h);

#endif