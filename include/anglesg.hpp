#ifndef _ANGLESG_
#define _ANGLESG_
/**
 * @file anglesg.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera anglesg
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\elements.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\global.hpp"
#include "..\include\IERS.hpp"
#include "..\include\geodetic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\timediff.hpp"
#include "..\include\hgibbs.hpp"
#include "..\include\gibbs.hpp"
#include "..\include\PrecMat.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\rpoly.hpp"
#include <tuple>
/**
 * @brief 
 * 
 * @param az1 
 * @param az2 
 * @param az3 
 * @param el1 
 * @param el2 
 * @param el3 
 * @param Mjd1 
 * @param Mjd2 
 * @param Mjd3 
 * @param Rs1 
 * @param Rs2 
 * @param Rs3 
 * @return tuple<Matrix&, Matrix&> 
 */
tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3);

#endif