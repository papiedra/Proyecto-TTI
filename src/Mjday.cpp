#include "..\include\Mjday.hpp"

double Mjday(int year, int month, int day, double hr=0.0, double min=0.0, double sec=0.0){
    double jd=367.0*year-floor((7*(year+floor((month+9)/12.0)))*0.25)+floor(275.0*month/9.0)+day+1721013.5+ (((sec / 60.0 + min) / 60.0 + hr) / 24.0);
    double Mjd=jd-2400000.5;
    return Mjd;
}