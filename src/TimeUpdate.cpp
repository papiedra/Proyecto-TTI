#include "../include/timeupdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt) {
    
    P = Phi * P * transpose(Phi) + Qdt;
    return P;
}