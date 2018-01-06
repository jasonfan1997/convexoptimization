#include<iostream>
#include<string>
#include <cmath>
#include <time.h>
#include <stdlib.h> 
#include "Matrix.h"
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono;

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
Matrix algorithm1withbaseline(Matrix& x,const double &iteration,const double &L,bool output,double baseline)
{
	if(output==1){
		for(int i=0;i<iteration&&F(x)>baseline;i++)
		{
				cout<<F(x)<<','<<endl;
				x=prox(x-(1/L)*grad_f(x),L);
		}
	}
	if(output==0){
		int i;
		for(i=0;i<iteration&&F(x)>baseline;i++)
		{
			x=prox(x-(1/L)*grad_f(x),L);
		}
		cout<<"Algorithm1: "<<i<<"steps"<<endl;
	}

	return x;
};
Matrix algorithm2withbaseline(Matrix& x,const double &iteration,const double &L,bool output,double baseline,double beta)
{
	if(output==1){
		Matrix y(N,1);
		Matrix temp(N,1);
		for(int i=0;i<iteration&&F(x)>baseline;i++)
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
		int i;
		for(i=0;i<iteration&&F(x)>baseline;i++)
		{
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			y=x+beta*(x-temp);
		}
		cout<<"Algorithm2: "<<i<<"steps"<<endl;
	}

	return x;
};
Matrix algorithm3withbaseline(Matrix& x,const double &iteration,const double &L,bool output,double baseline,double theta)
{
	if(output==1){
		Matrix y(N,1);
		y=x;
		Matrix temp(N,1);
		double beta;
		double temptheta;
		for(int i=0;i<iteration&&F(x)>baseline;i++)
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
		int i;
		for(i=0;i<iteration&&F(x)>baseline;i++)
		{
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			temptheta=theta;
			theta=(-theta+sqrt(5.0*(theta*theta)))/2;
			beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
			y=x+beta*(x-temp);
			
		}
		cout<<"Algorithm3: "<<i<<"steps"<<endl;
	}

	return x;
};
Matrix algorithm4withbaseline(Matrix& x,const double &iteration,const double &L,bool output,double baseline,double theta)
{
	if(output==1){
		Matrix y(N,1);
		Matrix z(N,1);
		z=x;
		for(int i=0;i<iteration&&F(x)>baseline;i++)
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
		int i;
		for(i=0;i<iteration&&F(x)>baseline;i++)
		{
			y=(1-theta)*x+theta*z;
			z=prox(z-(1/(theta*L))*grad_f(y),theta*L);
			x=prox(y-(1/L)*grad_f(y),L);
			theta=(-theta+sqrt(5.0*(theta*theta)))/2;
			
		}
		cout<<"Algorithm4: "<<i<<"steps"<<endl;
	}

	return x;
};
Matrix algorithm6withbaseline(Matrix& x,const double &iteration,const double &L,bool output,double baseline,double theta,double mu)
{
	if(output==1){
		Matrix y(N,1);
		y=x;
		Matrix temp(N,1);
		double beta;
		double temptheta;
		for(int i=0;i<iteration&&F(x)>baseline;i++)
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
		double temptheta;
		int i;
		for(i=0;i<iteration&&F(x)>baseline;i++)
		{
			temp=x;
			x=prox(y-(1/L)*grad_f(y),L);
			temptheta=theta;
			theta=(theta*theta)/(1-mu/L+theta*theta);
			beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
			y=x+beta*(x-temp);
			
		}
		cout<<"Algorithm6: "<<i<<"steps"<<endl;
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
	int iteration=200;
	double dummy;
	Matrix x(N,1);
	Matrix base(N,1);
	cout<<"Please enter L:";
	cin>>L;
	x=algorithm1(x,200,L,0);
	double baseline=F(x);
	cout<<"baseline="<<baseline<<" algorithm 1 with 200 steps"<<endl;
	x=base;
	int i;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	x=algorithm1withbaseline(x,iteration,L,0,baseline);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	cout<<"time used:"<<duration/1000000.0<<endl;
	x=base;
	t1 = high_resolution_clock::now();
	x=algorithm2withbaseline(x,iteration,L,0,baseline,0.5);
	t2 = high_resolution_clock::now();
	duration = duration_cast<microseconds>( t2 - t1 ).count();
	cout<<"time used:"<<duration/1000000.0<<endl;
	x=base;
	double theta=1;
	t1 = high_resolution_clock::now();
	x=algorithm3withbaseline(x,iteration,L,0,baseline,theta);
	t2 = high_resolution_clock::now();
	duration = duration_cast<microseconds>( t2 - t1 ).count();
	cout<<"time used:"<<duration/1000000.0<<endl;
	x=base;
	t1 = high_resolution_clock::now();
	x=algorithm4withbaseline(x,iteration,L,0,baseline,theta);
	t2 = high_resolution_clock::now();
	duration = duration_cast<microseconds>( t2 - t1 ).count();
	cout<<"time used:"<<duration/1000000.0<<endl;
	system("pause");
	
}
