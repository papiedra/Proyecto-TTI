/**
 * @file Position.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación Position
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\Position.hpp"
/**
 * @brief Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 * 
 * @param lon longitud
 * @param lat latitud
 * @param h altura
 * @return Matrix vector posición
 */
Matrix Position(double lon, double lat, double h){
    double R_equ = R_Earth;
    double f     = f_Earth;
    double e2=f*(2.0-f);
    double CosLat=cos(lat);
    double SinLat=sin(lat);
    double N=R_equ/sqrt(1.0-e2*SinLat*SinLat);
    Matrix r=zeros(3);
    r(1) = (N + h) * CosLat * cos(lon);
    r(2) = (N + h) * CosLat * sin(lon);
    r(3) = ((1.0 - e2) * N + h) * SinLat;
    return r;
}