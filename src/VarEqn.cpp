#include "..\include\vareqn.hpp"


Matrix& VarEqn(double x, Matrix yPhi) {

    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC,'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix   T, E, r, v, Phi, a, dfdy, Phip;
    Matrix P = PrecMat(MJD_J2000,AuxParam.Mjd_TT + x/86400);
    Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400);
    
    T = N * P;
    

    Matrix Pole=PoleMatrix(x_pole,y_pole) ;
    Matrix GHA=GHAMatrix(Mjd_UT1);
    E = Pole * GHA * T;
    r = extract_vector(transpose(yPhi),1,3);
    v = extract_vector(transpose(yPhi),4,6);
    Phi = zeros(6, 6);
    
    for (int j=1; j <= 6; j++) {
        assign_column(Phi, extract_vector(transpose(yPhi),6*j+1, 6*j+6),j);
    }
    
    a = transpose(AccelHarmonic ( transpose(r), E, AuxParam.n, AuxParam.m ));
    Matrix G = G_AccelHarmonic ( r, E, AuxParam.n, AuxParam.m );
    Matrix& yPhip = zeros(42,1);
    dfdy = zeros(6, 6);
    for (int i=1; i <= 3; i++) {
        for (int j=1; j <= 3; j++) {
            dfdy(i,j) = 0.0;                 
            dfdy(i+3,j) = G(i,j);            
            if ( i==j ) {
                dfdy(i,j+3) = 1;
            } else {
                dfdy(i,j+3) = 0;             
            }
            dfdy(i+3,j+3) = 0.0;             
        }
    }

    Phip = dfdy*Phi;

    for (int i=1; i <= 3; i++) {
        yPhip(i)   = v(i);                 
        yPhip(i+3) = a(i);                 
    }
    for (int i=1; i <= 6; i++) {
        for (int j=1; j<=6; j++) {
            yPhip(6*j+i) = Phip(i,j);     
        }
    }
    
    return yPhip;

}