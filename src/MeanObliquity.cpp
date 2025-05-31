/**
 * @file MeanObliquity.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación MeanObliquity
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\MeanObliquity.hpp"
#include "..\include\SAT_Const.hpp" 
/**
 * @brief Computes the mean obliquity of the ecliptic
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double Mean obliquity of the ecliptic [rad]
 */
double MeanObliquity(double Mjd_TT){
    double T=(Mjd_TT-MJD_J2000)/36525.0;
    double aux = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;
    return Rad * (aux / 3600.0);  
}

