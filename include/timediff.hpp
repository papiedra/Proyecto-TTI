#ifndef _TIMEDIFF_
#define _TIMEDIFF_
/**
 * @file timediff.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera timediff
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <cmath>
/**
 * @brief 
 * 
 * @param UT1_UTC 
 * @param TAI_UTC 
 * @param UT1_TAI 
 * @param UTC_GPS 
 * @param UT1_GPS 
 * @param TT_UTC 
 * @param GPS_UTC 
 */
void timediff(double UT1_UTC, double TAI_UTC,
              double& UT1_TAI, double& UTC_GPS,
              double& UT1_GPS, double& TT_UTC,
              double& GPS_UTC);

#endif