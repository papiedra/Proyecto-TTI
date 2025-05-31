#include "..\include\MeasUpdate.hpp"

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