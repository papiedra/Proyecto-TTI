#include "..\include\gmst.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
double gmst(double Mjd_UT1){
    double secs,MJD_J2000,gmstime,Mjd_0,UT1,T_0,T,gmst;
    secs = 86400.0;                       
    MJD_J2000 = 51544.5;

    Mjd_0 = floor(Mjd_UT1);
    UT1   = secs*(Mjd_UT1-Mjd_0);         
    T_0   = (Mjd_0  -MJD_J2000)/36525.0;
    T     = (Mjd_UT1-MJD_J2000)/36525.0;

    gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1+ (0.093104-6.2e-6*T)*T*T;    

    gmstime = 2*pi*(gmst/secs-std::floor(gmst/secs)); 
    return gmstime;
}