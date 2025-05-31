/**
 * @file unit.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación unit
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\unit.hpp"
/**
 * @brief this function calculates a unit vector given the original vector. if a zero vector is input, the vector is set to zero.
 * 
 * @param vec vector
 * @return Matrix& unit vector
 */
Matrix& unit(Matrix vec) {
    double small = 0.000001;
    double magv = norm(transpose(vec));
    Matrix* outvec = new Matrix(3);
    if ( magv > small ) {
        for (int i=1; i <= 3; i++) {
            (*outvec)(i)= vec(i)/magv;
        }
    } else {
        for (int i=1; i <= 3; i++) {
            (*outvec)(i)= 0.0;
        }
    }

    return *outvec;

}