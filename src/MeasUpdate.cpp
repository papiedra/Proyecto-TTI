#include "..\include\MeasUpdate.hpp"

std::tuple<Matrix, Matrix, Matrix> MeasUpdate(Matrix& x,Matrix& z,Matrix& g,Matrix& s,Matrix& G,Matrix& P,double n){
    int m=z.n_column;
    Matrix Inv_W=zeros(m,m);
    for (int i=1;i<=m;i++){
        Inv_W(i,i) = s(i)*s(i);    
    }
    Matrix K=P*transpose(G)*inv(Inv_W+(G*P*transpose(G)));
    Matrix x1=x+transpose(K*transpose(z-g));
    Matrix aux=z-g;
    P=(eye(n)-K*G)*P;
    return std::make_tuple(K,x1,P);
}