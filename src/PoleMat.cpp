#include "..\include\PoleMatrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include <cmath>

Matrix PoleMatrix(double xp, double yp){
    Matrix M=zeros(3,3);
    Matrix Rx=R_x(-yp);
    Matrix Ry=R_y(-xp);
    M=Ry*Rx;
    return M;  
}
