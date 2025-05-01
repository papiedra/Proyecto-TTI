#include "..\include\NutMatrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include <cmath>

Matrix NutMatrix(double Mjd_TT){
    double eps=MeanObliquity(Mjd_TT);
    std::tuple<double, double> aux=NutAngles(Mjd_TT);
    double deps=std::get<1>(aux);
    double dpsi=std::get<0>(aux);
    Matrix M=zeros(3,3);
    Matrix Rx1=R_x(-eps-deps);
    Matrix Rx2=R_x(eps);
    Matrix Rz=R_z(-dpsi);
    M=Rx1*Rz*Rx2;
    return M;  
}
