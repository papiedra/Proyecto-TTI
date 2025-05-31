/**
 * @file GHAMatrix.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación GHAMatrix
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\GHAMatrix.hpp"

/**
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system 
 * 
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return Matrix Greenwich Hour Angle matrix
 */
Matrix GHAMatrix(double Mjd_UT1){
   Matrix GHSmat=R_z(gast(Mjd_UT1));
   return GHSmat;
}
