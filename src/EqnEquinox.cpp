/**
 * @file EqnEquinox.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación EqnEquinox
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\EqnEquinox.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\MeanObliquity.hpp"
#include <cmath>
#include <tuple>
/**
 * @brief Computation of the equation of the equinoxes
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double Equation of the equinoxes
 */
double EqnEquinox(double Mjd_TT){
    double EqE;
    std::tuple<double, double> result = NutAngles(Mjd_TT);
    double dpsi=std::get<0>(result);
    EqE=dpsi*cos(MeanObliquity(Mjd_TT));
    return EqE;
}