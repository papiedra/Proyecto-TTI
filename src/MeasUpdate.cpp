/**
 * @file MeasUpdate.cpp
 * @author Pablo Piedrafita Sanromán
 * @brief Implementación MeasUpdate
 * @version 0.1
 * @date 2025-05-31
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "..\include\MeasUpdate.hpp"
/**
 * @brief Performs a measurement update step in a Kalman filte
 * 
 * @param x The current state vector 
 * @param z The measurement value
 * @param g The predicted measurement value
 * @param s The standard deviation of the measurement noise
 * @param G The measurement matrix
 * @param P The current state covariance matrix
 * @param n The dimension of the state vector
 * @return tuple<Matrix&, Matrix&, Matrix&> (The Kalman gain matrix, The updated state vector, The updated state covariance matrix)
 */
tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n) {
	
    double m = 1;
    Matrix Inv_W = zeros(m,m);
    Inv_W(1,1) = s*s;    
    Matrix& K = P*transpose(G)*inv(Inv_W+G*P*transpose(G));
	if(x.n_row==1){
		x=transpose(x);
	}
    Matrix& x2 = x + K*(z-g);
    Matrix& P2 = (eye(n)-K*G)*P;

    return tie(K, x2, P2);
}