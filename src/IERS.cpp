#include "..\include\matrix.hpp"
#include "..\include\IERS.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>
std::tuple<double,  double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp){
    


        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;


        if (interp =='l') {

            double mjd = (floor(Mjd_UTC));



            int i = 0;
            Matrix aux = extract_row(eop,4);
            for (int j = 1; j <= aux.n_column; j++) {
                if (aux(j) == mjd) {
                    i = j;
                    break;
                }
            }

    
            Matrix preeop = extract_column(eop,i);

    
            Matrix nexteop = extract_column(eop,i+1);
            
            double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
            double fixf = mfme/1440;

            x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
            y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
            UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
            LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
            dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
            deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
            dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
            dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
            TAI_UTC = preeop(13);
            
            x_pole  = x_pole/Arcs;  
            y_pole  = y_pole/Arcs;  
            dpsi    = dpsi/Arcs;
            deps    = deps/Arcs;
            dx_pole = dx_pole/Arcs; 
            dy_pole = dy_pole/Arcs; 
        } else if (interp =='n') {
            double mjd = (floor(Mjd_UTC));
            int i = 0;
            Matrix aux = extract_row(eop,4);
            for (int j = 1; j <= aux.n_column; j++) {
                if (aux(j) == mjd) {
                    i = j;
                    break;
                }
            }

            Matrix neweop = extract_column(eop,i);
            x_pole  = neweop(5)/Arcs;  
            y_pole  = neweop(6)/Arcs;  
            UT1_UTC = neweop(7);             
            LOD     = neweop(8);             
            dpsi    = neweop(9)/Arcs;
            deps    = neweop(10)/Arcs;
            dx_pole = neweop(11)/Arcs; 
            dy_pole = neweop(12)/Arcs; 
            TAI_UTC = neweop(13);           

        
    }
    return std::make_tuple(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}