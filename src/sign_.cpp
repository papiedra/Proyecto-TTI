/**
 * @file sign_.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación de sign_
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\sign_.hpp"
/**
 * @brief returns absolute value of a with sign of b
 * 
 * @param a número para calcular el abs
 * @param b número para el signo
 * @return double resultado
 */
double sign_(double a, double b){
    if(b>=0.0){
        return abs(a);
    }else{return -abs(a);}
}
