#ifndef _MJDAY_
#define _MJDAY_
/**
 * @file Mjday.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera Mjday
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
 * @param year 
 * @param month 
 * @param day 
 * @param hr 
 * @param min 
 * @param sec 
 * @return double 
 */
double Mjday(int year, int month, int day, double hr, double min, double sec);

#endif