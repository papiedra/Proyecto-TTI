#include "..\include\G_AccelHarmonic.hpp"


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
