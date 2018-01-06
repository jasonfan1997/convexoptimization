// The code needs you to finish it

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
using namespace std;
class Matrix
{

public:
	bool flag;
	int row,column;
	double **data;
	Matrix();
	Matrix(int r, int c);
	Matrix(const Matrix& m);
	Matrix& operator= (const Matrix& m);
	~Matrix();
	
	  // Member functions like
	double get_entry(int r, int c);
	int get_row_dimension() const; // the number of rows;
	int get_column_dimension() const;
	void transpose_prime();
		
    void row_switching(int i, int j);
    void row_multiplication(int i, double k);
    void row_addition(int i, int j, double k) ;

	void transpose();
  // Overloaded operators
    friend ostream& operator<<( ostream & cout, const  Matrix& t );
    friend istream & operator>> (istream& cin, Matrix& m);

 	friend Matrix operator+(const  Matrix& t,const  Matrix& u);
    friend Matrix operator-(const  Matrix& t,const  Matrix& u);
	friend Matrix operator*(const  Matrix& t,double a);
	friend Matrix operator*(double a,const  Matrix& t);
    friend Matrix operator*(const  Matrix& t,const  Matrix& u);

  // overload ">>", the input format:
  //   a list of double numbers delimited by ',', no space between numbers


	

};
#endif
