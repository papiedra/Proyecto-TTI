#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
	double **data;
	Matrix();
    // Parameterized constructor
    Matrix(const int n_row, const int n_column);
	Matrix(const int n);
	
	// Member operators
	double& operator () (const int row, const int column);
	double& operator () (const int n);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (Matrix &m);
	Matrix& operator / (Matrix &m);
	Matrix& operator = (Matrix &m);
	Matrix& operator + (double n);
	Matrix& operator - (double n);
	Matrix& operator * (double n);
	Matrix& operator / (double n);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);
Matrix& eye (int n);
Matrix& inv (Matrix &m);
Matrix& transpose (Matrix &m);
Matrix& zeros(const int n);
double norm(Matrix &m);
double dot(Matrix &m,Matrix &n);
Matrix& cross(Matrix &m,Matrix &n);
Matrix& extract_vector(Matrix &m, const int n1, const int n2);
Matrix& union_vector(Matrix &m, Matrix &n);
Matrix& extract_row(Matrix &m, const int n);
Matrix& extract_column(Matrix &m, const int n);
Matrix& assign_row(Matrix &m, Matrix &v, const int n);
Matrix& assign_column(Matrix &m, Matrix &v, const int n);
#endif