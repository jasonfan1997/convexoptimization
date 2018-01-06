// This program needs you to finish!
#include <iostream>
#include <cstring>
#include <string>
#include <cstdlib>
#include "Matrix.h"

using namespace std;

Matrix::Matrix(){
    flag=0;
    row=0;
    column=0;
  data=NULL;
}

Matrix::Matrix(int r, int c){
	if(r>0&&c>0){
    flag=0;
  row=r;
  column=c;
  data=new double *[r];
  for(int i=0;i<r;i++)
  {
      data[i]=new double[column];
    for(int j=0;j<c;j++)
    {
          data[i][j]=0;
    }
  }
}
else
{
	cout<<"Invalid number of rows or columns!"<<endl;
}
}
Matrix::Matrix(const Matrix& m)
{
    flag=m.flag;
	row=m.row;
	column=m.column;
	data=new double *[row];
  for(int i=0;i<row;i++)
  {
      data[i]=new double[column];
    for(int j=0;j<column;j++)
    {
          data[i][j]=0;
    }
  }
	for(int i=0;i<row;i++)
  {
    for(int j=0;j<column;j++)
    {
          data[i][j]=m.data[i][j];
    }
  }
}
void Matrix::row_switching(int i, int j)
{
	for(int k=0;k<column;k++)
	{
		double temp=data[i][k];
		data[i][k]=data[j][k];
		data[j][k]=temp;
	}
}
    void Matrix::row_multiplication(int i, double k)
	{
		for(int d=0;d<column;d++)
	{
		data[i][d]=k*data[i][d];

	}
	}
    void Matrix::row_addition(int i, int j, double k)
	{
		for(int d=0;d<column;d++)
	{
		data[i][d]+=k*data[j][d];

	}
	}
Matrix& Matrix::operator= (const Matrix& m){
row=m.row;
	column=m.column;
	flag=m.flag;
	data=new double *[row];
  for(int i=0;i<row;i++)
  {
     data[i]=new double[column];
    for(int j=0;j<column;j++)
    {
          data[i][j]=m.data[i][j];
    }
  }
  return *this;
}
Matrix::~Matrix(){
  // the destructor
  delete [] data;
}
   double Matrix::get_entry(int r, int c)
   {
   return data[r][c];
   }

int Matrix::get_row_dimension()const {
  // a const member function
  return row;
}

int Matrix::get_column_dimension() const
{
	return column;
}
 void Matrix::transpose()
 {
	Matrix m;
	m.row=column;
	m.column=row;
	m.data=new double *[m.row];
  for(int i=0;i<m.row;i++)
  {
      m.data[i]=new double[m.column];
    for(int j=0;j<m.column;j++)
    {
          m.data[i][j]=0;
    }
  }
	for(int i=0;i<row;i++)
  {
    for(int j=0;j<column;j++)
    {
          m.data[j][i]=data[i][j];
    }
  }
	row=m.row;
	column=m.column;
	for(int i=0;i<row;i++)
  {
    for(int j=0;j<column;j++)
    {
          data[i][j]=m.data[i][j];
    }
  }
 }
 void Matrix::transpose_prime()
 {
     	flag=1;
 }
ostream& operator<<( ostream & cout, const  Matrix& t )
 {
     if(t.flag==0){
	for (int i = 0; i < t.row; i++)
    {
		cout<<t.data[i][0];
         for (int j = 1; j < t.column; j++)
        {
            cout << " " <<t.data[i][j];
        }
        cout << endl;
    }
    return cout;
     }
 else{
    for (int i = 0; i < t.column; i++)
    {
		cout<<t.data[0][i];
         for (int j = 1; j < t.row; j++)
        {
            cout << " " <<t.data[j][i];
        }
        cout << endl;
    }
    return cout;
 }
 }

istream & operator>> (istream& cin, Matrix& m){
  // overload >> for input
   for (int i = 0; i < m.row; ++i)    {
       for (int j = 0; j < m.column; ++j)
        {
		cin>>m.data[i][j];
		}
		}
		return cin;
}
Matrix operator+(const  Matrix& t,const  Matrix& u)
{
	Matrix r;
	r.row=t.row;
	r.column=t.column;
	r.data=new double *[r.row];
  for(int i=0;i<r.row;i++)
  {
      r.data[i]=new double[r.column];
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=0;
    }
  }
		for(int i=0;i<r.row;i++)
  {
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=t.data[i][j]+u.data[i][j];
    }
  }
  return r;
}
Matrix operator-(const  Matrix& t,const  Matrix& u)
{
Matrix r;
	r.row=t.row;
	r.column=t.column;
	r.data=new double *[r.row];
  for(int i=0;i<r.row;i++)
  {
      r.data[i]=new double[r.column];
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=0;
    }
  }
		for(int i=0;i<r.row;i++)
  {
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=t.data[i][j]-u.data[i][j];
    }
  }
  return r;
}
Matrix operator*(double a,const  Matrix& t)
{
	Matrix r;
	r.row=t.row;
	r.column=t.column;
	r.data=new double *[r.row];
  for(int i=0;i<r.row;i++)
  {
      r.data[i]=new double[r.column];
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=0;
    }
  }
		for(int i=0;i<r.row;i++)
  {
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=a*t.data[i][j];
    }
  }
  return r;
}
Matrix operator*(const  Matrix& t,double a)
{
	return a*t;
	}
Matrix operator*(const  Matrix& t,const  Matrix& u)
{
	Matrix r;
	 if(t.flag==0&&u.flag==0)
  {
	r.row=t.row;
	r.column=u.column;
	r.data=new double *[r.row];
  for(int i=0;i<r.row;i++)
  {
      r.data[i]=new double[r.column];
    for(int j=0;j<r.column;j++)
    {
          r.data[i][j]=0;
    }
  }
	for(int i=0;i<r.row;i++)
  {
    for(int j=0;j<r.column;j++)
    {
         for(int k=0;k<t.column;k++)
  {
       r.data[i][j]+=t.data[i][k]*u.data[k][j];
  }
    }
  }
  }
   if(t.flag==1&&u.flag==1)
   {
       	Matrix c(t);
       	Matrix d(u);
       	c.flag=0;
       	d.flag=0;
       	r=d*c;
       	r.flag=1;
   }
   if(t.flag==0&&u.flag==1)
   {
       r.row=t.row;
       r.column=u.row;
       r.flag=0;
       	r.data=new double *[r.row];
  for(int i=0;i<r.row;i++)
  {
      r.data[i]=new double[r.column];
    for(int j=0;j<r.column;j++)
    {
           for(int k=0;k<t.column;k++)
  {
       r.data[i][j]+=t.data[i][k]*u.data[j][k];
  }
    }
  }
   }
     if(t.flag==1&&u.flag==0)
   {
       r.row=t.column;
       r.column=u.row;
       r.flag=0;
       	r.data=new double *[r.row];
  for(int i=0;i<r.row;i++)
  {
      r.data[i]=new double[r.column];
    for(int j=0;j<r.column;j++)
    {
           for(int k=0;k<r.column;k++)
  {
       r.data[i][j]+=t.data[k][i]*u.data[k][j];
  }
    }
  }
   }

return r;
}



