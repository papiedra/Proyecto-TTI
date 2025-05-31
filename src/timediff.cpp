/**
 * @file timediff.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación de timediff
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\timediff.hpp"

#include <cmath>
/**
 * @brief 
 * 
 * @param UT1_UTC Difference between UT1 and UTC time scales [s]
 * @param TAI_UTC Difference between TAI and UTC time scales [s]
 * @param UT1_TAI Difference between UT1 and TAI time scales [s] (out)
 * @param UTC_GPS Difference between UTC and GPS time scales [s] (out)
 * @param UT1_GPS Difference between UT1 and GPS time scales [s] (out)
 * @param TT_UTC Difference between TT and UTC time scales [s] (out)
 * @param GPS_UTC Difference between GPS and UTC time scales [s] (out)
 */
void timediff(double UT1_UTC, double TAI_UTC,
              double& UT1_TAI, double& UTC_GPS,
              double& UT1_GPS, double& TT_UTC,
              double& GPS_UTC) {

    const double TT_TAI  = 32.184;          
    const double GPS_TAI = -19.0;            
    const double TT_GPS  =  TT_TAI - GPS_TAI; 
    const double TAI_GPS = -GPS_TAI;        

    UT1_TAI = UT1_UTC - TAI_UTC;  
    double UTC_TAI = -TAI_UTC;   

    UTC_GPS = UTC_TAI - GPS_TAI;  
    UT1_GPS = UT1_TAI - GPS_TAI;  

    TT_UTC  = TT_TAI - UTC_TAI;   
    GPS_UTC = GPS_TAI - UTC_TAI;  
}