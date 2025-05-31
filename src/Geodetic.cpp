/**
 * @file Geodetic.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación Geodetic
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\Geodetic.hpp"
/**
 * @brief geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
 * 
 * @param r 
 * @return tuple<double, double, double> 
 */
tuple<double, double, double> Geodetic(Matrix r) {
    
    double R_equ = R_Earth;
    double f     = f_Earth;

    double epsRequ = 2.2204e-16*R_equ;        
    double e2      = f*(2.0-f);        
    double X = r(1);                   
    double Y = r(2);
    double Z = r(3);
    double rho2 = X*X + Y*Y;           
    double lon, lat, h;             
    if (norm(r)==0.0) {
        printf ( " invalid input in Geodetic constructor\n" );
        
        lon = 0.0;
        lat = 0.0;
        h   = -R_equ;
    }

    
    double dZ = e2*Z;
    double ZdZ,Nh,N;

    while(1) {
         ZdZ    =  Z + dZ;
         Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        double SinPhi =  ZdZ / Nh;                    
         N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        double dZ_new =  N*e2*SinPhi;
        if ( abs(dZ-dZ_new) < epsRequ ) {
            break;
        }
        dZ = dZ_new;
    }

    
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;
    return tie(lon, lat, h);

}