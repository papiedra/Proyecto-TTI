/**
 * @file G_AccelHarmonic.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación G_AccelHarmonic
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\G_AccelHarmonic.hpp"

/**
 * @brief Computes the gradient of the Earth's harmonic gravity field
 * 
 * @param r Satellite position vector in the true-of-date system
 * @param U Transformation matrix to body-fixed system
 * @param n_max Gravity model degree
 * @param m_max Gravity model order
 * @return Matrix Gradient (G=da/dr) in the true-of-date system
 */
Matrix G_AccelHarmonic(Matrix &r,Matrix &U,double n_max,double m_max){
   double d=1;
   Matrix G=zeros(3,3);
   Matrix dr=zeros(3);
   Matrix da(3,3);
   for(int i=1;i<=3;i++){
        dr=zeros(3);
        dr(i)=d;
        da = AccelHarmonic ( transpose(r+dr/2),U, n_max, m_max ) -  AccelHarmonic ( transpose(r-dr/2),U, n_max, m_max );
        assign_column(G, transpose(da), i);
   }
   return G;
}
