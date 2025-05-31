/**
 * @file Cheb3D.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación Cheb3D
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\Cheb3D.hpp"
/**
 * @brief Chebyshev approximation of 3-dimensional vectors
 * 
 * @param t 
 * @param N Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval
 * @param Cx Coefficients of Chebyshev polyomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polyomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polyomial (z-coordinate)
 * @return Matrix& 
 */
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz){
   
    if(t<Ta || Tb<t){
        cout<<"ERROR: Time out of range in Cheb3D"<<endl;
        exit(EXIT_FAILURE);
    }
    
    double tau=(2*t-Ta-Tb)/(Tb-Ta);
    Matrix f1=zeros(3);
    Matrix f2=zeros(3);
    Matrix old_f1=zeros(3);
    for(int i=N;i>=2;i--){
        old_f1=f1;
        Matrix aux(3);
        aux(1) = Cx(i); aux(2) = Cy(i); aux(3) = Cz(i);
        f1=f1*(2*tau)-f2+aux;
        f2=old_f1;
    }
    Matrix aux1(3);
    aux1(1) = Cx(1); aux1(2) = Cy(1); aux1(3) = Cz(1);
    Matrix &ChebApp=((f1*tau)-f2)+aux1;
    return ChebApp;
}