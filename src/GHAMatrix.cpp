#include "..\include\GHAMatrix.hpp"


Matrix GHAMatrix(double Mjd_UT1){
   Matrix GHSmat=R_z(gast(Mjd_UT1));
   return GHSmat;
}
