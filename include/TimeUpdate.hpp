#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_
/**
 * @file TimeUpdate.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabecera TimeUpdate
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "matrix.hpp"
/**
 * @brief 
 * 
 * @param P 
 * @param Phi 
 * @param Qdt 
 * @return Matrix& 
 */
Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt=0.0);



#endif