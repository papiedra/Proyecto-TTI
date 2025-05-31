/**
 * @file gast.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación gast
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\gast.hpp"
#include "..\include\gmst.hpp"
#include "..\include\EqnEquinox.hpp"
/**
 * @brief Greenwich Apparent Sidereal Time
 * 
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return double GAST in [rad]
 */
double gast(double Mjd_UT1){
    double gstime = fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * M_PI);
    return gstime;
}
