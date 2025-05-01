#include "../include/AccelHarmonic.hpp"
#include "../include/global.hpp"
#include "../include/SAT_Const.hpp"
#include "../include/Legendre.hpp"

Matrix& AccelHarmonic(Matrix& r, Matrix E, int n_max,int m_max){
    double r_ref, gm, d, latgc, lon, dUdr, dUdlatgc, dUdlon, q3, q2, q1, b1, b2, b3, r2xy, ax, ay, az;
    Matrix r_bf, a_bf;
    a_bf=zeros(3);
    r_ref=6378.1363e3;
    gm=398600.4415e9;
    r_bf=E*r;
    r_bf=transpose(r_bf);
    d=norm(r_bf);
    latgc=asin(r_bf(3)/d);
    lon=atan2(r_bf(2),r_bf(1));
    auto [pnm,dpnm]=legendre(n_max,m_max,latgc);
    dUdr = 0;
    dUdlatgc = 0;
    dUdlon = 0;
    q3 = 0; q2 = q3; q1 = q2;
    for (int n=0; n <= n_max; n++) {
        b1 = (-gm/pow(d, 2))*pow((r_ref/d),n)*(n+1);
        b2 =  (gm/d)*pow(r_ref/d,n);
        b3 =  (gm/d)*pow(r_ref/d,n);
        for (int m=0; m <= m_max; m++) {
            q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
        }
        dUdr = dUdr+ q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;
        q3 = 0; q2 = q3; q1 = q2;
        
    }
    r2xy = pow(r_bf(1),2)+pow(r_bf(2),2);
    
    ax = (1/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
    ay = (1/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
    az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/pow(d,2)*dUdlatgc;
    a_bf(1)=ax; a_bf(2)=ay; a_bf(3)=az;
    a_bf = transpose(a_bf);
    Matrix& aux = transpose(E)*a_bf;
    return aux;
}