#include "..\include\matrix.hpp"

Matrix::Matrix(){
	this->n_row=0;
	this->n_column=0;
	this->data=nullptr;
}
Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

Matrix::Matrix(const int n) {
    if (n <= 0) {
		cout << "Matrix create: error number\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = 1;
	this->n_column = n;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	this->data[0] = (double *) malloc(n_column*sizeof(double));
	
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

double& Matrix::operator () (const int n) {
	if (n <= 0 || n > this->n_column*this->n_row) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n-1)/this->n_column][(n-1)%this->n_column];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator * (Matrix &m){
	if (this->n_column != m.n_row) {
		cout << "Matrix mult: El número de columnas de la primera matriz debe ser igual al número de filas de la segunda\n";
        exit(EXIT_FAILURE);
	}
	Matrix *m_aux = new Matrix(this->n_row, m.n_column);
	for(int i=1;i<=this->n_row;i++){
		for(int j=1;j<=m.n_column;j++){
			(*m_aux)(i,j)=0;
			for(int k=1;k<=this->n_column;k++){
				(*m_aux)(i,j)+=(*this)(i,k)*m(k,j);
			}
		}
	}
	return *m_aux;
}

Matrix& Matrix::operator / (Matrix &m){
	if (this->n_column != m.n_row) {
		cout << "Matrix div: tamaños incompatibles para división\n";
        exit(EXIT_FAILURE);
	}
	Matrix invertida=inv(m);
	return *this * invertida;
}

Matrix& Matrix::operator = (Matrix &m)
{
	this->n_row = m.n_row;
	this->n_column = m.n_column;

	this->data = (double **) malloc(m.n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix assignment: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < m.n_row; i++) {
		this->data[i] = (double *) malloc(m.n_column*sizeof(double));
		for (int j = 0; j < this->n_column; j++) {
			this->data[i][j]=m.data[i][j];
		}
	}
	
	return *this;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& Matrix::operator + (double n)
{
	Matrix *aux=new Matrix(this->n_row,this->n_column);
    for (int i = 1; i <= n_row; i++)
        for (int j = 1; j <= n_column; j++){
            (*aux)(i,j) = (*this)(i,j)+n;
		}
    return *aux;
}

Matrix& Matrix::operator - (double n)
{
	Matrix *aux=new Matrix(this->n_row,this->n_column);
    for (int i = 1; i <= n_row; i++)
        for (int j = 1; j <= n_column; j++)
		(*aux)(i,j) = (*this)(i,j)-n;
 
    return *aux;
}

Matrix& Matrix::operator * (double n)
{
	Matrix *aux=new Matrix(this->n_row,this->n_column);
    for (int i = 1; i <= n_row; i++)
        for (int j = 1; j <= n_column; j++)
		(*aux)(i,j) = (*this)(i,j)*n;
 
    return *aux;
}

Matrix& Matrix::operator / (double n)
{
	if(abs(n)<1e-10){
		cout << "Matrix double div: División por 0\n";
        exit(EXIT_FAILURE);
	}
	Matrix *aux=new Matrix(this->n_row,this->n_column);
    for (int i = 1; i <= n_row; i++)
        for (int j = 1; j <= n_column; j++)
		(*aux)(i,j) = (*this)(i,j)/n;
    return *aux;
}

Matrix& zeros(const int n_row, const int n_column) {
	
	Matrix *m_aux = new Matrix(n_row, n_column);
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	return (*m_aux);
}

Matrix& eye(int n) {
	if (n<=0) {
		cout << "Matrix eye: error debe ser de al menos tamaño 1\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(n, n);
	
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
			if(i==j){
				(*m_aux)(i,j) = 1;
			}else{
				(*m_aux)(i,j) =0;
			}
		}
	}
	
	return *m_aux;
}

Matrix& inv(Matrix &m) {
	if (m.n_column!=m.n_row) {
		cout << "Matrix inv: error debe ser cuadrada\n";
        exit(EXIT_FAILURE);
	}
	Matrix *inv=new Matrix(m);
	Matrix *aux=&eye(m.n_column);
	for(int i=1;i<=m.n_column;i++){
		double pivot=(*inv)(i,i);
		if(fabs(pivot)<1e-10){
			std::cout << "Matrix inv: matriz no invertible\n";
            exit(EXIT_FAILURE);
		}
		for(int j=1;j<=m.n_row;j++){
			(*inv)(i,j)/=pivot;
			(*aux)(i,j)/=pivot;
		}
		for(int k=1;k<=m.n_row;k++){
			if(k!=i){
				double f=(*inv)(k,i);
				for(int j=1;j<=m.n_row;j++){
					(*inv)(k,j)-=(*inv)(i,j)*f;
					(*aux)(k,j)-=(*aux)(i,j)*f;
				}
			}
		}
	}
	return *aux;
	
}
Matrix& transpose(Matrix &m) {
	Matrix *aux=new Matrix(m.n_column,m.n_row);
	for(int i=1;i<=m.n_column;i++){
		for(int j=1;j<=m.n_row;j++){
			(*aux)(j,i)=m(i,j);
		}
	}
	return *aux;
	
}

Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n);
	
	for(int i = 1; i <= n; i++) {
		(*m_aux)(n)=0;
	}
	
	return (*m_aux);
}

double norm(Matrix &m) {
	double res=0;
	for(int i=1;i<=m.n_column;i++){
		res+=pow(m(i),2);
	}
	double raiz=sqrt(res);
	return raiz;
}
double dot(Matrix &m,Matrix &n) {
	double res=0;
	for(int i=1;i<=m.n_column;i++){
		res+=m(i)*n(i);
	}
	
	return res;
}
Matrix& cross(Matrix &m,Matrix &n) {
	if (m.n_column!=3 || n.n_column!=3 || 1!=m.n_row || 1!=n.n_row ) {
		cout << "Matrix cross: error deben ser vectores de R3\n";
        exit(EXIT_FAILURE);
	}
	Matrix *aux = new Matrix(3);
	(*aux)(1)=m(2)*n(3)-m(3)*n(2);
	(*aux)(2)=m(3)*n(1)-m(1)*n(3);
	(*aux)(3)=m(1)*n(2)-m(2)*n(1);
	
	return (*aux);
}

Matrix& extract_vector(Matrix &m,const int n1, const int n2) {
	if (n1>n2 || n1<=0 || n2<=0 || n2>m.n_column) {
		cout << "Extract vector: Index out of bounds\n";
        exit(EXIT_FAILURE);
	}
	int tam=n2-n1+1;
	Matrix *aux = new Matrix(tam);
	for(int j=1;j<=tam;j++){
		(*aux)(j)=m(n1+j-1);
	}
	return (*aux);
}

Matrix& union_vector(Matrix &m,Matrix &n) {
	if (m.n_row>1 || n.n_row>1) {
		cout << "Union vector: Los parámetros no son vectores\n";
        exit(EXIT_FAILURE);
	}
	int tam=m.n_column+n.n_column;
	Matrix *aux = new Matrix(tam);
	for(int j=1;j<=m.n_column;j++){
		(*aux)(j)=m(j);
	}
	for(int i=m.n_column+1;i<=tam;i++){
		(*aux)(i)=n(i-m.n_column);
	}
	return (*aux);
}

Matrix& extract_row(Matrix &m,const int n) {
	if (n>m.n_row) {
		cout << "Extract row: Index out of bounds\n";
        exit(EXIT_FAILURE);
	}
	int tam=m.n_column;
	Matrix *aux = new Matrix(tam);
	for(int j=1;j<=tam;j++){
		(*aux)(j)=m(n,j);
	}
	return (*aux);
}

Matrix& extract_column(Matrix &m,const int n) {
	if (n>m.n_column) {
		cout << "Extract column: Index out of bounds\n";
        exit(EXIT_FAILURE);
	}
	int tam=m.n_row;
	Matrix *aux = new Matrix(tam);
	for(int j=1;j<=tam;j++){
		(*aux)(j)=m(j,n);
	}
	return (*aux);
}

Matrix& assign_row(Matrix &m, Matrix &v, const int n) {
	if (n>m.n_row || v.n_column!=m.n_column || v.n_row!=1) {
		cout << "Assing row: Index out of bounds\n";
        exit(EXIT_FAILURE);
	}
	Matrix *aux=&m;
	for(int j=1;j<=m.n_column;j++){
		(*aux)(n,j)=v(j);
	}
	return (*aux);
}

Matrix& assign_column(Matrix &m, Matrix &v, const int n ) {
	if (n>m.n_column || v.n_column!=m.n_row || v.n_row!=1) {
		cout << "Assing column: Index out of bounds\n";
        exit(EXIT_FAILURE);
	}
	Matrix *aux=&m;
	for(int j=1;j<=m.n_row;j++){
		(*aux)(j,n)=v(j);
	}
	return (*aux);
}