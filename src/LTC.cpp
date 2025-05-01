#include "..\include\LTC.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include <cmath>

Matrix LTC(double lon, double lat){
    Matrix M=zeros(3,3);
    Matrix Ry=R_y(-1.0*lat);
    Matrix Rz=R_z(lon);
    M=Ry*Rz;
    double Aux;
    for (int j=1;j<=3;j++){
        Aux=M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;

    }
    return M;
    
}
