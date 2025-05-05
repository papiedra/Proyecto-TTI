#include "..\include\gast.hpp"
#include "..\include\gmst.hpp"
#include "..\include\EqnEquinox.hpp"

double gast(double Mjd_UT1){
    double gstime = fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * M_PI);
    return gstime;
}
