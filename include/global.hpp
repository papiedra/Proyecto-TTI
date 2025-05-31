#ifndef _GLOBAL_
#define _GLOBAL_
/**
 * @file global.hpp
 * @author Pablo Piedrafita Sanromán
 * @brief Cabeceras de los métodos que cargan ficheros y declaraciones de las variables extern
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "matrix.hpp"
#include "SAT_Const.hpp"
#include "..\include\Mjday.hpp"
#include <cmath>
#include <string.h>

typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern Param AuxParam;
extern Matrix Obs;
/**
 * @brief 
 * 
 * @param c 
 */
void eop19620101(int c);
/**
 * @brief 
 * 
 */
void GGM03S();
/**
 * @brief 
 * 
 */
void DE430Coeff();
/**
 * @brief 
 * 
 */
void auxparam();
/**
 * @brief 
 * 
 * @param f 
 */
void GEOS3(int f);
#endif