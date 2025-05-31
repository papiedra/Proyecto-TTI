/**
 * @file EFK_GEOS3.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Metodo main
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Position.hpp"
#include "..\include\accel.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\LTC.hpp"
#include "..\include\anglesg.hpp"
#include "..\include\gmst.hpp"
#include <iostream>

using namespace std;
/**
 * @brief Initial Orbit Determination using Gauss and Extended Kalman Filter methods
 * 
 * @return int 
 */
int main() {
	eop19620101(21413); 
    GGM03S();
    DE430Coeff();
    auxparam();
    int nObs = 46;
    GEOS3(nObs);

    double  Mjd_TT, Mjd_UT1, theta, Dist, Azim, Elev;
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC, t, t_old;
    Matrix  Y_old,  dAdY, dEdY, dDds, dDdY, K, Y0, Y_true;

    double sigma_range = 92.5;          
    double sigma_az = 0.0224*Rad; 
    double sigma_el = 0.0139*Rad; 

   
    double lat = Rad*21.5748;     
    double lon = Rad*(-158.2706); 
    double alt = 300.20;           
    Matrix aux= Position(lon, lat, alt);
    Matrix Rs = transpose(aux);

    double Mjd1 = Obs(1,1);
    double Mjd2 = Obs(9,1);
    double Mjd3 = Obs(18,1);

    //auto [r2,v2] = anglesg(Obs(1,2),Obs(9,2),Obs(18,2),Obs(1,3),Obs(9,3),Obs(18,3),
                  //Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

     Matrix r2(3, 1);
     r2(1,1)=6221397.62857869;
     r2(2,1)=2867713.77965738;
     r2(3,1)=3006155.98509949;
     Matrix v2(3,1);
     v2(1,1)=4645.04725161806;
     v2(2,1)=-2752.21591588204;
     v2(3,1)=-7507.99940987031;
    // Comentario auxiliar para cuando no tenía implementada la iteración extra.
     
    Matrix Y0_apr = transpose(union_vector(transpose(r2), transpose(v2)));
    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = Obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;

    int n_eqn  = 6;
    
    Matrix Y = transpose(DEInteg(Accel,0,-(Obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr));   
    Matrix P = zeros(6, 6);
    for (int i=1; i <= 3; i++) {
        P(i,i)=1e8;
    }
    for (int i=4; i <= 6; i++) {
        P(i,i)=1e3;
    }
    Matrix LT = LTC(lon,lat);

    Matrix yPhi = zeros(42,1);
    Matrix Phi  = zeros(6, 6);

      
    t = 0;
  
    for (int i=1; i <= nObs; i++) {    

        t_old = t;
        Y_old = Y;
        
        Mjd_UTC = Obs(i,1);                       
        t       = (Mjd_UTC-Mjd0)*86400.0;         

        tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
		double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
 
        for (int ii=1; ii <= 6; ii++) {
            yPhi(ii) = Y_old(ii);
            for (int j=1; j <= 6; j++) {  
                if (ii==j) {
                    yPhi(6*j+ii) = 1; 
                } else {
                    yPhi(6*j+ii) = 0;
                }
            }
        }
        if (yPhi.n_row==1) {
            yPhi = transpose(yPhi);
        }

        yPhi = transpose(DEInteg(VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi));

        for (int j=1; j <= 6; j++) {
            assign_column(Phi,extract_vector(yPhi,6*j+1, 6*j+6),j);
        }
        if (Y_old.n_row==1) {
            Y_old=transpose(Y_old);
        }
        Y = transpose(DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old));
       
        theta = gmst(Mjd_UT1);                    
        Matrix U = R_z(theta);
        Matrix r = transpose(extract_vector(Y,1, 3)); 
        Matrix s = LT*(U*r-Rs);                          
        P = TimeUpdate(P, Phi);
        std::tuple<double, double, Matrix, Matrix> result = AzElPa(s);
        Azim=std::get<0>(result);
        Elev=std::get<1>(result);
        Matrix dAds=std::get<2>(result);
        Matrix dEds=std::get<3>(result);     

        dAdY = union_vector(dAds*LT*U, zeros(1,3));
        if (Y.n_row==1) {
            Y=transpose(Y);
        }
        tie(K, Y, P) = MeasUpdate ( Y, Obs(i,2), Azim, sigma_az, dAdY, P, 6 );

    
        r = transpose(extract_vector(transpose(Y),1, 3)); 
        s = LT*(U*r-Rs);                            
        std::tuple<double, double, Matrix, Matrix> result1 = AzElPa(s);
        Azim=std::get<0>(result1);
        Elev=std::get<1>(result1);
        Matrix dAds1=std::get<2>(result1);
        Matrix dEds1=std::get<3>(result1);       
        dEdY = union_vector(dEds1*LT*U, zeros(1,3));
    
        tie(K, Y, P) = MeasUpdate ( Y, Obs(i,3), Elev, sigma_el, dEdY, P, 6 );
     
        r = transpose(extract_vector(Y,1, 3)); 
        s = LT*(U*r-Rs);                          
        Dist = norm(s); 
		dDds = transpose(s/Dist);        

        dDdY = union_vector(dDds*LT*U,zeros(1,3));
    
        tie(K, Y, P) = MeasUpdate ( Y, Obs(i,4), Dist, sigma_range, dDdY, P, 6 );
    }
	
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Obs(46,1),'l');
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;
    Y0 = DEInteg (Accel,0,-(Obs(46,1)-Obs(1,1))*86400.0,1e-13,1e-6,6,Y);
    
    Y_true = zeros(6);
    Y_true(1) = 5753.173e3;
    Y_true(2) = 2673.361e3;
    Y_true(3) = 3440.304e3;
    Y_true(4) = 4.324207e3;
    Y_true(5) = -1.924299e3;
    Y_true(6) = -5.728216e3;
    
    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n", Y0(1)-Y_true(1));
    printf("dY%10.1f [m]\n", Y0(2)-Y_true(2)); 
    printf("dZ%10.1f [m]\n", Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n", Y0(4)-Y_true(4));
    printf("dVy%8.1f [m/s]\n", Y0(5)-Y_true(5));
    printf("dVz%8.1f [m/s]\n", Y0(6)-Y_true(6));

}