/**
 * @file AccelPointMass.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación AccelPointMass
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\AccelPointMass.hpp"
/**
 * @brief Computes the perturbational acceleration due to a point mass
 * 
 * @param r Satellite position vector 
 * @param s Point mass position vector
 * @param GM Gravitational coefficient of point mass
 * @return Matrix& Acceleration (a=d^2r/dt^2)
 */
Matrix& AccelPointMass(Matrix &r, Matrix &s, double GM){
    Matrix &d=r-s;
    Matrix &a=((d / pow(norm(d), 3)) + (s / pow(norm(s), 3)))*(-GM);
    return a;
}