/**
 * @file TimeUpdate.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación TimeUpdate
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "../include/timeupdate.hpp"
/**
 * @brief Performs the time update step for a covariance matrix in a Kalman filter
 * 
 * @param P Reference to the current covariance matrix to be updated
 * @param Phi Reference to the state transition matrix
 * @param Qdt Process noise covariance scaled by the time step
 * @return Matrix& Matriz de covarianza después de la actualización de tiempo
 */
Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt) {
    
    P = Phi * P * transpose(Phi) + Qdt;
    return P;
}