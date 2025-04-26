#include "..\include\matrix.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_mult_01() {
    int f1 = 4;
    int c1 = 3;
	int f2 = 3;
    int c2 = 5;
	
	Matrix A(f1, c1);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3; 
	A(2,1) = 4; A(2,2) = 5; A(2,3) = 6; 
	A(3,1) = 7; A(3,2) = 8; A(3,3) = 9; 
	A(4,1) = 0; A(4,2) = 1; A(4,3) = 2; 
	
	Matrix B(f2, c2);
	B(1,1) = 0; B(1,2) = 1; B(1,3) = 2; B(1,4) = 3; B(1,5) = 4;
	B(2,1) = 5; B(2,2) = 6; B(2,3) = 7; B(2,4) = 8; B(2,5) = 9;
	B(3,1) = 0; B(3,2) = 1; B(3,3) = 2; B(3,4) = 3; B(3,5) = 4;

	Matrix C(f1, c2);
	C(1,1) = 10; C(1,2) = 16; C(1,3) = 22; C(1,4) = 28; C(1,5) = 34;
	C(2,1) = 25; C(2,2) = 40; C(2,3) = 55; C(2,4) = 70; C(2,5) = 85;
	C(3,1) = 40; C(3,2) = 64; C(3,3) = 88; C(3,4) = 112; C(3,5) = 136;
	C(4,1) = 5; C(4,2) = 8; C(4,3) = 11; C(4,4) = 14; C(4,5) = 17;
    Matrix res=A*B;
    _assert(m_equals(C, res, 1e-10));
    
    return 0;
}
int m_div_01() {
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	Matrix B(3, 3);
	B(1,1) = 1; B(1,2) = 0; B(1,3) = 0;  
	B(2,1) = 0; B(2,2) = 1; B(2,3) = 0; 
	B(3,1) = 0; B(3,2) = 0; B(3,3) = 1; 
 
	Matrix res(3,3);
	res(1,1) = 1; res(1,2) = 0; res(1,3) = 2;  
	res(2,1) = 0; res(2,2) = 1; res(2,3) = 0; 
	res(3,1) = 0; res(3,2) = 0; res(3,3) = 1; 
    _assert(m_equals(A/B, res, 1e-10));
    
    return 0;
}

int m_eye_01() {
	
	Matrix A(4, 4);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0; 
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; A(2,4) = 0; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; A(3,4) = 0; 
	A(4,1) = 0; A(4,2) = 0; A(4,3) = 0; A(4,4) = 1; 
	
    Matrix res=eye(4);
    _assert(m_equals(A, res, 1e-10));
    
    return 0;
}
int m_inv_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
 
	Matrix res(3,3);
	res(1,1) = 1; res(1,2) = 0; res(1,3) = -2;  
	res(2,1) = 0; res(2,2) = 1; res(2,3) = 0; 
	res(3,1) = 0; res(3,2) = 0; res(3,3) = 1; 
    _assert(m_equals(inv(A), res, 1e-10));
    
    return 0;
}
int m_eq_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	Matrix res=A;
    _assert(m_equals(A, res, 1e-10));
    
    return 0;
}
int m_trans_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 3; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	Matrix res=transpose(A);
	Matrix B(3, 3);
	B(1,1) = 1; B(1,2) = 0; B(1,3) = 0;  
	B(2,1) = 3; B(2,2) = 1; B(2,3) = 0; 
	B(3,1) = 2; B(3,2) = 3; B(3,3) = 1; 
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_nsum_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 3; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	double n=3;
	Matrix res=A+3;
	Matrix B(3, 3);
	B(1,1) = 4; B(1,2) = 6; B(1,3) = 5;  
	B(2,1) = 3; B(2,2) = 4; B(2,3) = 6; 
	B(3,1) = 3; B(3,2) = 3; B(3,3) = 4; 
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_nsub_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 3; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	double n=3;
	Matrix res=A-3;
	Matrix B(3, 3);
	B(1,1) = -2; B(1,2) = 0; B(1,3) = -1;  
	B(2,1) = -3; B(2,2) = -2; B(2,3) = 0; 
	B(3,1) = -3; B(3,2) = -3; B(3,3) = -2; 
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_nmul_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 3; A(1,3) = 2;  
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1; 
	double n=3;
	Matrix res=A*3;
	Matrix B(3, 3);
	B(1,1) = 3; B(1,2) = 9; B(1,3) = 6;  
	B(2,1) = 0; B(2,2) = 3; B(2,3) = 9; 
	B(3,1) = 0; B(3,2) = 0; B(3,3) = 3; 
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_ndiv_01() {
	
	Matrix A(3, 3);
	A(1,1) = 6; A(1,2) = 3; A(1,3) = 9;  
	A(2,1) = 0; A(2,2) = 6; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 6; 
	double n=3;
	Matrix res=A/n;
	Matrix B(3, 3);
	B(1,1) = 2; B(1,2) = 1; B(1,3) = 3;  
	B(2,1) = 0; B(2,2) = 2; B(2,3) = 1; 
	B(3,1) = 0; B(3,2) = 0; B(3,3) = 2; 
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_extn_01() {
	
	Matrix A(3, 3);
	A(1,1) = 6; A(1,2) = 3; A(1,3) = 9;  
	A(2,1) = 0; A(2,2) = 6; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 6; 
	double n=A(6);
	double res=3;
    _assert(m_equals(n, res, 1e-10));
    
    return 0;
}
int m_nzeros_01() {
    int c = 4;
	
	Matrix A(c);
	A(1) = 0; A(2) = 0; A(3) = 0; A(4) = 0;
	
	Matrix B = zeros(4);
	
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}
int m_norm_01() {
    int c = 4;
	
	Matrix A(c);
	A(1) = 1; A(2) = 2; A(3) = 3; A(4) = 4;
	
	double res=norm(A);
	double B = sqrt(30);
    
    _assert(m_equals(res, B, 1e-10));
    
    return 0;
}
int m_dot_01() {
    int c = 4;
	
	Matrix A(c);
	A(1) = 1; A(2) = 2; A(3) = 3; A(4) = 4;
	Matrix B(c);
	B(1) = 4; B(2) = 3; B(3) = 2; B(4) = 1;
	
	double res=20;
	double aux = dot(A,B);
    
    _assert(m_equals(res, aux, 1e-10));
    
    return 0;
}
int m_cross_01() {
    int c = 3;
	
	Matrix A(c);
	A(1) = 1; A(2) = 2; A(3) = 3; 
	Matrix B(c);
	B(1) = 3; B(2) = 2; B(3) = 1; 
	
	Matrix C(c);
	C(1) = -4; C(2) = 8; C(3) = -4; 
	Matrix res = cross(A,B);
    
    _assert(m_equals(res, C, 1e-10));
    
    return 0;
}
int m_extract_01() {
	
	Matrix A(6);
	A(1) = 1; A(2) = 2; A(3) = 3; A(4) = 1; A(5) = 2; A(6) = 3; 
	Matrix B=extract_vector(A,2,4);
	
	Matrix C(3);
	C(1) = 2; C(2) = 3; C(3) = 1; 
    
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}
int m_union_01() {
	int c = 3;
	Matrix A(c);
	A(1) = 1; A(2) = 2; A(3) = 3; 
	Matrix B(c);
	B(1) = 4; B(2) = 5; B(3) = 6; 
	
	Matrix C(6);
	C(1) = 1; C(2) = 2; C(3) = 3; C(4) = 4; C(5) = 5; C(6) = 6;
    Matrix res = union_vector(A,B);
    _assert(m_equals(res, C, 1e-10));
    
    return 0;
}
int m_extractrow_01() {
	
	Matrix A(3, 3);
	A(1,1) = 6; A(1,2) = 3; A(1,3) = 9;  
	A(2,1) = 0; A(2,2) = 6; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 6;
	
	Matrix C(3);
	C(1) = 0; C(2) = 6; C(3) = 3; 
    Matrix res=extract_row(A,2);
    _assert(m_equals(res, C, 1e-10));
    
    return 0;
}
int m_extractcol_01() {
	
	Matrix A(3, 3);
	A(1,1) = 6; A(1,2) = 3; A(1,3) = 9;  
	A(2,1) = 0; A(2,2) = 6; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 6;
	
	Matrix C(3);
	C(1) = 3; C(2) = 6; C(3) = 0; 
    Matrix res=extract_column(A,2);
    _assert(m_equals(res, C, 1e-10));
    
    return 0;
}
int m_assrow_01() {
	
	Matrix A(3, 3);
	A(1,1) = 6; A(1,2) = 3; A(1,3) = 9;  
	A(2,1) = 0; A(2,2) = 6; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 6; 
	Matrix C(3);
	C(1) = 3; C(2) = 6; C(3) = 0;
	Matrix res=assign_row(A,C,3);
	Matrix B(3, 3);
	B(1,1) = 6; B(1,2) = 3; B(1,3) = 9;  
	B(2,1) = 0; B(2,2) = 6; B(2,3) = 3; 
	B(3,1) = 3; B(3,2) = 6; B(3,3) = 0;  
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_asscol_01() {
	
	Matrix A(3, 3);
	A(1,1) = 6; A(1,2) = 3; A(1,3) = 9;  
	A(2,1) = 0; A(2,2) = 6; A(2,3) = 3; 
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 6; 
	Matrix C(3);
	C(1) = 3; C(2) = 6; C(3) = 0;
	Matrix res=assign_column(A,C,3);
	Matrix B(3, 3);
	B(1,1) = 6; B(1,2) = 3; B(1,3) = 3;  
	B(2,1) = 0; B(2,2) = 6; B(2,3) = 6; 
	B(3,1) = 0; B(3,2) = 0; B(3,3) = 0;  
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_AccelPointMass_01() {
	Matrix r(3);
	r(1) = 1; r(2) = 2; r(3) = 3; 
	Matrix s(3);
	s(1) = 3; s(2) = 2; s(3) = 1; 
	double GM=10;
	Matrix B=AccelPointMass(r,s,GM);
	Matrix res(3);
	res(1)=0.3112 ; res(2)= -0.3818; res(3)= -1.0748;
    _assert(m_equals(B, res, 1e-4)); // Es el nivel de precisión que me da matlab
    
    return 0;
}
int m_Cheb3D_01() {
	int N=5;
	double Ta=0.0; double Tb=10.0; double t=5.0;
	Matrix Cx(5);
	Cx(1) = 1; Cx(2) = 2; Cx(3) = 3; Cx(4) = 4; Cx(5) = 5; 
	Matrix Cy(5);
	Cy(1) = 1.5; Cy(2) = 2.5; Cy(3) = 3.5; Cy(4) = 4.5; Cy(5) = 5.5; 
	Matrix Cz(5);
	Cz(1) = 2; Cz(2) = 3; Cz(3) = 4; Cz(4) = 5; Cz(5) = 6; 
	Matrix B=Cheb3D(t,N,Ta,Tb,Cx,Cy,Cz);
	Matrix res(3);
	res(1)=3 ; res(2)= 3.5; res(3)= 4;
    _assert(m_equals(B, res, 1e-10));
    
    return 0;
}
int m_EccAnom_01() {
	double M=1.5;
	double e=0.5;
	double res=EccAnom(M,e);
	double B=1.9622;
    _assert(abs(res-B<1e-10));
    
    return 0;
}
int m_Frac_01() {
	double x=1.5;
	double res=Frac(x);
	double B=1;
    _assert(abs(res-B<1e-10));
    
    return 0;
}
int m_MeanObliquity_01() {
	double Mjd_TT=58000;
	double res=Frac(Mjd_TT);
	double B=0.000241073;
    _assert(abs(res-B<1e-10));
    
    return 0;
}
int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_zeros_01);
	_verify(m_mult_01);
	_verify(m_div_01);
	_verify(m_eye_01);
	_verify(m_inv_01);
	_verify(m_eq_01);
	_verify(m_trans_01);
	_verify(m_nsum_01);
	_verify(m_nsub_01);
	_verify(m_nmul_01);
	_verify(m_ndiv_01);
	_verify(m_extn_01);
	_verify(m_nzeros_01);
	_verify(m_norm_01);
	_verify(m_dot_01);
	_verify(m_cross_01);
	_verify(m_extract_01);
	_verify(m_union_01);
	_verify(m_extractrow_01);
	_verify(m_extractcol_01);
	_verify(m_assrow_01);
	_verify(m_asscol_01);
	_verify(m_AccelPointMass_01);
	_verify(m_Cheb3D_01);
	_verify(m_EccAnom_01);
	_verify(m_Frac_01);
	_verify(m_MeanObliquity_01);
    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
