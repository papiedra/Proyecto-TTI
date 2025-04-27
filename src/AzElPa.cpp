#include "..\include\matrix.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>
std::tuple<double, double, Matrix, Matrix> AzElPa(Matrix& s){
    double rho = sqrt(s(1) * s(1) + s(2) * s(2));

    double Az = atan2(s(1), s(2)); 

    if (Az < 0.0) {
        Az += pi2;
    }

    double El = atan(s(3) / rho); 

    // Partials
    Matrix dAds=zeros(3);
    dAds(1) = s(2) / (rho * rho);
    dAds(2) = -s(1) / (rho * rho);
    dAds(3) = 0.0;
    Matrix dEds=zeros(3);
    dEds(1) = -s(1) * s(3) / rho;
    dEds(2) = -s(2) * s(3) / rho;
    dEds(3) = rho;

    // Normalize dEds by the dot product of s
    double dot_product = dot(s,s);
    Matrix res = dEds/dot_product;
    return std::make_tuple(Az,El,dAds,dEds);
}