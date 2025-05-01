
#include "..\include\EqnEquinox.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\MeanObliquity.hpp"
#include <cmath>
#include <tuple>

double EqnEquinox(double Mjd_TT){
    double EqE;
    std::tuple<double, double> result = NutAngles(Mjd_TT);
    double dpsi=std::get<0>(result);
    EqE=dpsi*cos(MeanObliquity(Mjd_TT));
    return EqE;
}