#include<iostream>
#include<string>
#include <cmath>
#include <time.h>
#include<vector>
#include <stdlib.h> 
#include "Matrix.h"
#include <random>
#include <fstream>
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
			theta=(-theta+sqrt(pow(theta,4)+4*pow(theta,2)))/2;
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
			theta=(-theta+sqrt(pow(theta,4)+4*pow(theta,2)))/2;
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

void APG1(Matrix &x,Matrix &y,double &theta,double &beta,double temptheta,Matrix &temp,double &L)
{
	temp=x;
	x=prox(y-(1/L)*grad_f(y),L);
	temptheta=theta;
	theta=(-theta+sqrt(pow(theta,4)+4*pow(theta,2)))/2;
	beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
	y=x+beta*(x-temp);
}
void APG2(Matrix &x,Matrix &y,Matrix &z,double &theta,double &L)
{
	y=(1-theta)*x+theta*z;
	z=prox(z-(1/(theta*L))*grad_f(y),theta*L);
	x=prox(y-(1/L)*grad_f(y),L);
	theta=(-theta+sqrt(pow(theta,4)+4*pow(theta,2)))/2;
}
void APG3(Matrix &x,Matrix &y,Matrix &z,double &A,double &a,Matrix &alpha,Matrix &x0,double &L)
{
	a=(2+sqrt(4+8*A*L))/(2*L);
	y=(A/(A+a))*x+(a/(A+a))*z;
	x=prox(y-(1/L)*grad_f(y),L);
	alpha=alpha+a*grad_f(x);
	for(int i=0;i<x.get_row_dimension();i++)
	{
		if((x0.data[i][0]-alpha.data[i][0])>A)
		{
			z.data[i][0]=x0.data[i][0]-alpha.data[i][0]-A;
		}
		else if(abs(x0.data[i][0]-alpha.data[i][0])<A)
		{
			z.data[i][0]=0;
		}
		else if((x0.data[i][0]-alpha.data[i][0])<-A)
		{
			z.data[i][0]=x0.data[i][0]-alpha.data[i][0]+A;
		}
	}
	
}


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
	cout<<"5.APG3"<<endl;
	cout<<"6.MAPG1"<<endl;
	cin>>selection;
	cout<<"Please enter the number of iteration:"<<endl;
	cin>>iteration;
	cout<<"Please enter the L:"<<endl;
	cin>>L;
//	cout<<"Output?(0 or 1)"<<endl;
//	cin>>output;
	output=1;
	if(selection==1){
		x=algorithm1(x,iteration,L,output);
	}
	if(selection==2){
		x=algorithm2(x,iteration,L,output);
	}
	if(selection==3)
	{
		ofstream myfile;
      	myfile.open ("APG1.csv");
		cout<<"Please type in theta: "<<endl;
		double theta;
		cin>>theta;
		double beta=0;
		Matrix y(N,1);
		Matrix temp(N,1);
		double temptheta=theta;
		y=x;
		temp=x;
		for(int i=0;i<iteration;i++)
		{
			APG1(x,y,theta,beta,temptheta,temp,L);
			myfile <<i+1<<","<<F(x)<<","<<F(y)<<endl;
		}	
		myfile.close();
		
	}
	if(selection==4)
	{
		ofstream myfile;
      	myfile.open ("APG2.csv");
		cout<<"Please type in theta: "<<endl;
		double theta;
		cin>>theta;
		Matrix y(N,1);
		Matrix z(N,1);
		y=x;
		z=x;
		for(int i=0;i<iteration;i++)
		{
			APG2(x,y,z,theta,L);
			myfile<<i+1<<","<<F(x)<<","<<F(y)<<endl;
		}	
		myfile.close();
	}
	if(selection==5)
	{
		ofstream myfile;
      	myfile.open ("APG3.csv");
		Matrix y(N,1);
		Matrix x0(N,1);
		Matrix z(N,1);
		y=x;
		z=x;
		x0=x;
		double a=0;
		double A=0;
		Matrix alpha;
		alpha=x;
		for(int i=0;i<iteration;i++)
		{
			APG3(x,y,z,A,a,alpha,x0,L);
			myfile<<i+1<<","<<F(x)<<","<<F(y)<<endl;
		}
		myfile.close();
	}
	if(selection==6)
	{
		x=algorithm6(x,iteration,L,output);
	}
	system("pause");
	
}
