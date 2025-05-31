/**
 * @file R_y.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación R_y
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\R_y.hpp"
/**
 * @brief Genera la matriz de rotación para el ángulo pasado como parámetro en torno al eje y
 * 
 * @param angle angle of rotation [rad]
 * @return Matrix result
 */
Matrix R_y(double angle){
    double C=cos(angle);
    double S=sin(angle);
    Matrix rotmat=zeros(3,3);
    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;
    return rotmat;
}