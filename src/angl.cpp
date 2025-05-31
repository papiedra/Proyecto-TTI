/**
 * @file angl.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación angl
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\angl.hpp"
/**
 * @brief angle between the two vectors
 * 
 * @param vec1 vector 1
 * @param vec2 vector 2
 * @return double angle between the two vectors, pi to pi
 */
double angl(Matrix vec1, Matrix vec2){
    double small     = 0.00000001;
    double undefined = 999999.1;
    double magv1 = norm(vec1);
    double magv2 = norm(vec2);
    double theta  = 0.0;
	double temp;
    if (magv1*magv2 > pow(small,2)) {
        temp= dot(vec1,vec2) / (magv1*magv2);
        if (abs( temp ) > 1.0) {
            double sign = temp > 0 ? 1 : -1;
            temp= sign* 1.0;
        }
         theta= acos( temp );
    } else {
        theta= undefined;
    }

    return theta;
}