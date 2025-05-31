/**
 * @file LTC.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación LTC
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\LTC.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include <cmath>
/**
 * @brief  Transformation from Greenwich meridian system to local tangent coordinates
 * 
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @return Matrix Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 */
Matrix LTC(double lon, double lat){
    Matrix M=zeros(3,3);
    Matrix Ry=R_y(-1.0*lat);
    Matrix Rz=R_z(lon);
    M=Ry*Rz;
    double Aux;
    for (int j=1;j<=3;j++){
        Aux=M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;

    }
    return M;
    
}
