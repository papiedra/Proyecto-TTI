/**
 * @file PoleMat.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación PoleMatrix
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\PoleMatrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include <cmath>
/**
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 * 
 * @param xp coordenada x
 * @param yp coordenada y
 * @return Matrix Pole matrix
 */
Matrix PoleMatrix(double xp, double yp){
    Matrix M=zeros(3,3);
    Matrix Rx=R_x(-yp);
    Matrix Ry=R_y(-xp);
    M=Ry*Rx;
    return M;  
}
