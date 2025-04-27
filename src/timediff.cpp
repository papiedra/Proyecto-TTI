#include "..\include\timediff.hpp"

#include <cmath>

std::tuple<double,double,double,double,double> timediff(double UT1_UTC, double TAI_UTC){

    const double TT_TAI  = 32.184;          
    const double GPS_TAI = -19.0;            
    const double TT_GPS  =  TT_TAI - GPS_TAI; 
    const double TAI_GPS = -GPS_TAI;

    double UT1_TAI = UT1_UTC - TAI_UTC;  
    double UTC_TAI = -TAI_UTC;    

    double UTC_GPS = UTC_TAI - GPS_TAI;  
    double UT1_GPS = UT1_TAI - GPS_TAI;  

    double TT_UTC  = TT_TAI - UTC_TAI;   
    double GPS_UTC = GPS_TAI - UTC_TAI;  
    return std::make_tuple(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
}