// NOTA: Algunos valores para los test los he obtenido de chatgp o pidiéndolos a otros compañeros con el fin de obtener valores decentes de test.
#include "..\include\matrix.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\SAT_Const.hpp" 
#include "..\include\MeanObliquity.hpp" 
#include "..\include\Mjday_TBD.hpp"
#include "..\include\Position.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\global.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
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
int AccelPointMass_01() {
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
int Cheb3D_01() {
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
int EccAnom_01() {
	double M=1.5;
	double e=0.5;
	double res=EccAnom(M,e);
	double B=1.96219;
    _assert(abs(res-B)<1e-6);
    
    return 0;
}
int Frac_01() {
	double x=1.5;
	double res=Frac(x);
	double B=0.5;
    _assert(abs(res-B)<1e-6);
    
    return 0;
}
int MeanObliquity_01() {
	double Mjd_TT=58000;
	double res=MeanObliquity(Mjd_TT);
	double B=0.409053;
    _assert(abs(res-B)<1e-6);
    
    return 0;
}
int Mjday_01() {
	int yr=2025;
	int mon=4;
	int day=27;
	double hr=12;
	double min=0;
	double sec=0;
	double res=Mjday(yr,mon,day,hr,min, sec);
	double B=60792.5;
    _assert(abs(res-B)<1e-6);
    
    return 0;
}
int Mjday_TBD_01() {
	double Mjd_TT=51544.5;
	double res=Mjday_TBD(Mjd_TT);
	double B=51544.4999999988;
    _assert(abs(res-B)<1e-6);
    
    return 0;
}
int Position_01() {
    Matrix Pos = Position(0,1,1);
	Matrix A(3);  
    A(1) = 3454319.32320978910000000000;
    A(2) = 0;
    A(3) = 5343769.28735907750000000000;
	
	_assert(m_equals(A, Pos, 1e-10));

    return 0;
}
int R_x_01() {
	Matrix A = R_x(1.0);
	
    Matrix res(3, 3);  
    res(1,1) = 1.0000    ; res(1,2) = 0; res(1,3) =  0;
    res(2,1) = 0; res(2,2) = 0.54030230586814; res(2,3) =  0.841470984807897;
    res(3,1) = 0; res(3,2) = -0.841470984807897 ; res(3,3) = 0.54030230586814;

	
    _assert(m_equals(A, res, 1e-10));  

    return 0;
}
int R_y_01() {
	Matrix A = R_y(1.0);
	
    Matrix res(3, 3);  
    res(1,1) =0.54030230586814     ; res(1,2) = 0; res(1,3) =  -0.841470984807897;
    res(2,1) = 0; res(2,2) = 1; res(2,3) =  0;
    res(3,1) = 0.841470984807897; res(3,2) = 0 ; res(3,3) = 0.54030230586814;

	
    _assert(m_equals(A, res, 1e-10));  

    return 0;
}
int R_z_01() {
	Matrix A = R_z(1.0);
	
    Matrix res(3, 3);  
    res(1,1) =0.54030230586814     ; res(1,2) = 0.841470984807897 ; res(1,3) =  0;
    res(2,1) = -0.841470984807897 ; res(2,2) = 0.540302305868141; res(2,3) =  0;
    res(3,1) = 0; res(3,2) = 0 ; res(3,3) = 1;

	
    _assert(m_equals(A, res, 1e-10));  

    return 0;
}
int sign__01() {
	double a=5.0;
	double b=-2;
	double res=sign_(a,b);
	double aux=-5.0;
	
    _assert(res==aux);  

    return 0;
}
int timediff_01() {
	double UT1_UTC = -0.3341;
    double TAI_UTC = 37;
	auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC]=timediff(UT1_UTC,TAI_UTC);

    _assert(fabs(UT1_TAI + 37.3341) < 1e-10);
    _assert(fabs(UTC_GPS + 18.0) < 1e-10);
    _assert(fabs(UT1_GPS + 18.3341) < 1e-10);
    _assert(fabs(TT_UTC - 69.184) < 1e-10);
    _assert(fabs(GPS_UTC - 18.0) < 1e-10);

    return 0;
}
int AzElPa_01() {
	
    Matrix r(3);
	r(1)=1; r(2)=2; r(3)=3;
	double expected = 0.463647609000806; 
	double expected2 = 0.930274014115472;
	std::tuple<double, double, Matrix, Matrix> result = AzElPa(r);

	_assert(fabs(std::get<0>(result) - expected) < 1e-10);
	_assert(fabs(std::get<1>(result) - expected2) < 1e-10);
	return 0;
}
int IERS_01() {
	std::tuple<double, double, double, double, double, double, double, double, double> result = IERS(eopdata, 37670,'n');
	_assert(fabs(std::get<0>(result) - (-1.338037278494208e-07)) < 1e-10);
	_assert(fabs(std::get<1>(result) - 1.058353113998928e-06) < 1e-10);
	_assert(fabs(std::get<2>(result) - 0.030535300000000) < 1e-10);

	return 0;
}
int Legendre_01() {
	std::tuple<Matrix, Matrix> result = legendre(1, 2, 1.0);
	Matrix P = std::get<0>(result);
	Matrix dP = std::get<1>(result);
	Matrix expected(2,3);
	expected(1,1) = 1.0;    expected(1,2) = 0.0;       expected(1,3) = 0.0;
	expected(2,1) = 1.457470498782296; expected(2,2) = 0.935831045210238;    expected(2,3) = 0.0;
	_assert(m_equals(P, expected, 1e-10));
	return 0;
}
int NutAngles_01() {
	double Mjd_TT = 6.067604166666651e04;
	double expected = 9.723503682287755e-07;
	double expected2 = 4.120518762261807e-05;
	std::tuple<double, double> result = NutAngles(Mjd_TT);

	_assert(fabs(std::get<0>(result) - expected) < 1e-10);
	_assert(fabs(std::get<1>(result) - expected2) < 1e-10);
	return 0;
}
int TimeUpdate_01() {
	Matrix P(2, 2);
    P(1, 1) = 1;
    P(1, 2) = 2;
    P(2, 1) = 3;
    P(2, 2) = 4;

    Matrix Phi(2, 2);
    Phi(1, 1) = 2;
    Phi(1, 2) = 3;
    Phi(2, 1) = 1;
    Phi(2, 2) = 2;

    double Qdt = 1;

    Matrix P_updated = TimeUpdate(P, Phi, Qdt);
	Matrix Res(2, 2);
    Res(1, 1) = 71;
    Res(1, 2) = 44;
    Res(2, 1) = 45;
    Res(2, 2) = 28;
	
	_assert(m_equals(Res, P_updated, 1e-10));
	
	return 0;
}
int all_tests()
{
	
	eop19620101(10);
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
	_verify(AccelPointMass_01);
	_verify(Cheb3D_01);
	_verify(EccAnom_01);
	_verify(Frac_01);
	_verify(MeanObliquity_01);
	_verify(Mjday_01);
	_verify(Mjday_TBD_01);
	_verify(Position_01);
	_verify(R_x_01);
	_verify(R_y_01);
	_verify(R_z_01);
	_verify(sign__01);
	_verify(timediff_01);
	_verify(AzElPa_01);
	_verify(IERS_01);
	_verify(Legendre_01);
	_verify(NutAngles_01);
	_verify(TimeUpdate_01);
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
