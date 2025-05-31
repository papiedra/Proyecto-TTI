/**
 * @file NutMtrix.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación NutMatrix
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\NutMatrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include <cmath>
/**
 * @brief Transformation from mean to true equator and equinox
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Matrix Nutation matrix
 */
Matrix NutMatrix(double Mjd_TT){
    double eps=MeanObliquity(Mjd_TT);
    std::tuple<double, double> aux=NutAngles(Mjd_TT);
    double deps=std::get<1>(aux);
    double dpsi=std::get<0>(aux);
    Matrix M=zeros(3,3);
    Matrix Rx1=R_x(-eps-deps);
    Matrix Rx2=R_x(eps);
    Matrix Rz=R_z(-dpsi);
    M=Rx1*Rz*Rx2;
    return M;  
}
