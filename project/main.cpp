#include<iostream>
#include<string>
#include <cmath>
#include <time.h>
#include <stdlib.h> 
#include "Matrix.h"
#include <random>
#include <algorithm>
using namespace std;


int N=1000; //define the problem
int M=1000;
double L; //learing rate
Matrix A(N,M);  
Matrix b(N,1);

double squarevector(const Matrix &M) //taking square of the matrix
{
	double result=0;
	for(int i=0;i<M.get_row_dimension();i++)
	{
		result+=M.data[i][0]*M.data[i][0];
	}
	return result;
};
double innerproduct(const Matrix &M,const Matrix &N) //taking inner product of the matrix
{
	double result=0;
	for(int i=0;i<M.get_row_dimension();i++)
	{
		result+=M.data[i][0]*N.data[i][0];
	}
	return result;
};

double f(Matrix &x)
{
	return 0.5*squarevector(A*x-b);
};


double g(const Matrix &x)
{
	double result=0;
	for(int i=0;i<x.get_row_dimension();i++)
	{
		result+=abs(x.data[i][0]);
	}
	return result;
};
double F(Matrix &x)
{
	return f(x)+g(x);
};

Matrix grad_f(Matrix &x)
{
	return A*x-b;
};
Matrix prox(const Matrix& x,const double &L)
{
	Matrix result=x;
	for(int i=0;i<x.row;i++)
	{
		if(x.data[i][0]>=0)
		result.data[i][0]=max(abs(x.data[i][0])-1/L,0.0);
		else 
		result.data[i][0]=-max(abs(x.data[i][0])-1/L,0.0);
	}
	return result;
};
Matrix algorithm1(Matrix& x,const double &iteration,const double &L,bool output)
{
	if(output==1){
		for(int i=0;i<iteration;i++)
		{
				cout<<F(x)<<','<<endl;
				x=prox(x-(1/L)*grad_f(x),L);
		}
	}
	if(output==0){
		for(int i=0;i<iteration;i++)
		{
			x=prox(x-(1/L)*grad_f(x),L);
		}
	}

	return x;
};
Matrix algorithm2(Matrix& x,const double &iteration,const double &L,bool output)
{
	if(output==1){
		Matrix y(N,1);
		Matrix temp(N,1);
		double beta;
		cout<<"Please choose a beta between 0 and 1"<<endl;
		cin>>beta;
		for(int i=0;i<iteration;i++)
		{
			cout<<F(x)<<','<<endl;
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			y=x+beta*(x-temp);
		}
	}
	if(output==0){
		Matrix y(N,1);
		Matrix temp(N,1);
		double beta;
		cout<<"Please choose a beta between 0 and 1"<<endl;
		cin>>beta;
		for(int i=0;i<iteration;i++)
		{
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			y=x+beta*(x-temp);
		}
	}

	return x;
};
Matrix algorithm3(Matrix& x,const double &iteration,const double &L,bool output)
{
	if(output==1){
		Matrix y(N,1);
		y=x;
		Matrix temp(N,1);
		double beta;
		double theta;
		double temptheta;
		cout<<"Please choose a theta between 0 and 1"<<endl;
		cin>>theta;
		for(int i=0;i<iteration;i++)
		{
			cout<<F(x)<<','<<endl;
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			temptheta=theta;
			theta=(-theta+sqrt(5.0*(theta*theta)))/2;
			beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
			y=x+beta*(x-temp);
			
		}
	}
	if(output==0){
		Matrix y(N,1);
		y=x;
		Matrix temp(N,1);
		double beta;
		double theta;
		double temptheta;
		cout<<"Please choose a theta between 0 and 1"<<endl;
		cin>>theta;
		for(int i=0;i<iteration;i++)
		{
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			temptheta=theta;
			theta=(-theta+sqrt(5.0*(theta*theta)))/2;
			beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
			y=x+beta*(x-temp);
			
		}
	}

	return x;
};
Matrix algorithm4(Matrix& x,const double &iteration,const double &L,bool output)
{
	if(output==1){
		Matrix y(N,1);
		Matrix z(N,1);
		z=x;
		double theta;
		cout<<"Please choose a theta between 0 and 1"<<endl;
		cin>>theta;
		for(int i=0;i<iteration;i++)
		{
			cout<<F(x)<<','<<endl;
			y=(1-theta)*x+theta*z;
			z=prox(z-(1/(theta*L))*grad_f(y),theta*L);
			x=prox(y-(1/L)*grad_f(y),L);
			theta=(-theta+sqrt(5.0*(theta*theta)))/2;
			
		}
	}
	if(output==0){
		Matrix y(N,1);
		Matrix z(N,1);
		z=x;
		double theta;
		cout<<"Please choose a theta between 0 and 1"<<endl;
		cin>>theta;
		for(int i=0;i<iteration;i++)
		{
			y=(1-theta)*x+theta*z;
			z=prox(z-(1/(theta*L))*grad_f(y),theta*L);
			x=prox(y-(1/L)*grad_f(y),L);
			theta=(-theta+sqrt(5.0*(theta*theta)))/2;
			
		}
	}

	return x;
};
Matrix algorithm6(Matrix& x,const double &iteration,const double &L,bool output)
{
	if(output==1){
		Matrix y(N,1);
		y=x;
		Matrix temp(N,1);
		double beta;
		double theta;
		double mu;
		double temptheta;
		cout<<"Please choose a theta between 0 and 1"<<endl;
		cin>>theta;
		cout<<"Please enter the value of mu"<<endl;
		cin>>mu;
		for(int i=0;i<iteration;i++)
		{
			cout<<F(x)<<','<<endl;
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			temptheta=theta;
			theta=(theta*theta)/(1-mu/L+theta*theta);
			beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
			y=x+beta*(x-temp);
			
		}
	}
	if(output==0){
		Matrix y(N,1);
		y=x;
		Matrix temp(N,1);
		double beta;
		double theta;
		double mu;
		double temptheta;
		cout<<"Please choose a theta between 0 and 1"<<endl;
		cin>>theta;
		cout<<"Please enter the value of mu"<<endl;
		cin>>mu;
		for(int i=0;i<iteration;i++)
		{
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			temptheta=theta;
			theta=(theta*theta)/(1-mu/L+theta*theta);
			beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
			y=x+beta*(x-temp);
			
		}
	}

	return x;
};

int main()
{
	srand (time(NULL)); //set random seed
	for(int i=0;i<A.row;i++)  //init the matrix
	{
		
	b.data[i][0]=rand() % 100/10.0;
		for(int j=0;j<A.column;j++)
		{
			if(i==j)
				A.data[i][j]=rand() % 100/10.0;
		}
	}	
	int selection;
	int iteration;
	bool output;
	Matrix x(N,1);
	cout<<"Please choose the alogrithm:"<<endl;
	cout<<"1.Proximal gardient"<<endl;
	cout<<"2.Extraploted proximal gardient"<<endl;
	cout<<"3.APG1"<<endl;
	cout<<"4.APG2"<<endl;
	cout<<"6.MAPG1"<<endl;
	cin>>selection;
	cout<<"Please enter the number of iteration:"<<endl;
	cin>>iteration;
	cout<<"Please enter the L:"<<endl;
	cin>>L;
	cout<<"Output?(0 or 1)"<<endl;
	cin>>output;
	if(selection==1){
		x=algorithm1(x,iteration,L,output);
	}
	if(selection==2){
		x=algorithm2(x,iteration,L,output);
	}
	if(selection==3)
	{
		x=algorithm3(x,iteration,L,output);
	}
	if(selection==4)
	{
		x=algorithm4(x,iteration,L,output);
	}
	if(selection==6)
	{
		x=algorithm6(x,iteration,L,output);
	}
	system("pause");
	
}
