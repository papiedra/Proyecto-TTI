/**
 * @file Frac.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación Frac
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\Frac.hpp"
/**
 * @brief Fractional part of a number (y=x-[x])
 * 
 * @param x 
 * @return double 
 */
double Frac(double x){
    return x-floor(x);
}
