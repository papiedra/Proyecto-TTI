#include "..\include\EccAnom.hpp"
#include "..\include\SAT_Const.hpp" 
double EccAnom(double M, double e){
   int maxit=15;
   int i=1;
   M=fmod(M, 2.0 * M_PI);
   double E;
   if(e<0.8){
    E=M;
   }else{
    E=M_PI;
   }
   double f=E-e*sin(E)-M;
   E=E-f/(1-e*cos(E));
   while(fabs(f)> eps){
        f=E-e*sin(E)-M;
        E = E - f / ( 1.0 - e*cos(E) );
         i = i+1;
        if (i==maxit){
        cout<<" convergence problems in EccAnom"<<endl;
        exit(EXIT_FAILURE);
    }
   }
   return E;
}