/**
 * @file Mjday.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación Mjday
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\Mjday.hpp"
/**
 * @brief Calculates the modified julian date given the date and the time
 * 
 * @param year 
 * @param month 
 * @param day 
 * @param hr universal time hour  
 * @param min universal time min 
 * @param sec universal time sec 
 * @return double Modified julian date
 */
double Mjday(int year, int month, int day, double hr=0.0, double min=0.0, double sec=0.0){
    double jd=367.0*year-floor((7*(year+floor((month+9)/12.0)))*0.25)+floor(275.0*month/9.0)+day+1721013.5+ (((sec / 60.0 + min) / 60.0 + hr) / 24.0);
    double Mjd=jd-2400000.5;
    return Mjd;
}