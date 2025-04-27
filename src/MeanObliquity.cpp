#include "..\include\MeanObliquity.hpp"
#include "..\include\SAT_Const.hpp" 

double MeanObliquity(double Mjd_TT){
    double T=(Mjd_TT-MJD_J2000)/36525.0;
    double aux = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;
    return Rad * (aux / 3600.0);  
}

