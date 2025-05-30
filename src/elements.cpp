#include "..\include\elements.hpp"
tuple<double, double, double, double, double, double, double> elements(Matrix y) {
        
    double pi2 = pi2;
    
    Matrix r =extract_vector(y,1,3);                                        
    Matrix v =extract_vector(y,4,6);                                        
    
    Matrix h = cross(r,v);                                    
    double magh = norm(h);
    double p = magh*magh/GM_Earth;
    double H = norm(h);
    
    double Omega = atan2 ( h(1), -h(2) );                     
    Omega = fmod(Omega,pi2);
    if (Omega < 0) {
        Omega += pi2;
    }
    double i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) );         
    double u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );       
    
    double R  = norm(r);                                               
    
    double a = 1.0/(2.0/R-dot(v,v)/GM_Earth);               
    
    double eCosE = 1.0-R/a;                                            
    double eSinE = dot(r,v)/sqrt(GM_Earth*a);                 
    
    double e2 = eCosE*eCosE +eSinE*eSinE;
    double e  = sqrt(e2);                                       
    double E  = atan2(eSinE,eCosE);                              
    
    double M  = fmod(E-eSinE,pi2);                            
    if (M < 0) {
        M+=pi2;
    }
    
    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);           
    
    double omega = fmod(u-nu,pi2);                            
    if (omega < 0) {
        omega += pi2;
    }
    
    return tie(p, a, e, i, Omega, omega, M);
    
}