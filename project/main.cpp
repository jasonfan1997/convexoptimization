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

int N= 1000; //define the problem
int M= 100;
Matrix A(N,M);
Matrix b(N,1);
int Max= 100000000;
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
Matrix PG(Matrix& x,const int &iteration,const double &L)
{
	Matrix f_value(iteration,1);
	for (int i=0; i< iteration; i++){
		x=prox(x-(1/L)*grad_f(x),L);
		f_value.data[i][0]= f(x);
	}
	return f_value;	
};
void apg1(Matrix &x,Matrix &y,double &theta,double &mu,const double &L)
{
	Matrix temp;
	temp =x;
	x=prox(y-(1/L)*grad_f(y),L);
	double temptheta=theta;
	theta=(-pow(theta,2)+ mu/L+sqrt(pow(pow(theta,2)- mu/L,2)+4*pow(theta,2)))/2;
	double beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
	y=x+beta*(x-temp);
}
void apg2(Matrix &x,Matrix &z,double &theta,double &mu, const double &L, bool choice)
{
	Matrix y;
	y=(1-theta)*x+theta*z;
	z=prox(z-(1/(theta*L))*grad_f(y),theta*L);
	if (choice==0){
		x=prox(y-(1/L)*grad_f(y),L);
	}
	else{
		x= (1- theta)*x+ theta*z;
	}
	theta=(-pow(theta,2)+ mu/L+sqrt(pow(pow(theta,2)- mu/L,2)+4*pow(theta,2)))/2;
}
void apg3(Matrix &x,Matrix &z,Matrix &y,double &A,Matrix &alpha,Matrix &x0,double &mu,const double &L)
{
	double a=(2*(1+ mu*A)/L+sqrt(4*pow((1+ mu*A)/L,2)+8*(1+ mu*A)*A/L))/2;
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
	A= A+ a;
}
Matrix APG1(Matrix& x,const int &iteration,double &mu, const double &L)
{
	Matrix f_value(iteration,2);
	Matrix y;
	y= x;
	double theta= 1.0;
	for (int i=0; i< iteration; i++){
		apg1(x,y,theta,mu,L);
		f_value.data[i][0]= F(x);
		f_value.data[i][1]= F(y);
	}
	return f_value;
};
Matrix APG2(Matrix& x,const int &iteration,double &mu, const double &L)
{
	Matrix f_value(iteration,2);
	Matrix z;
	z= x;
	double theta= 1.0;
	for (int i=0; i< iteration; i++){
		apg2(x,z,theta,mu,L,0);
		f_value.data[i][0]= F(x);
		f_value.data[i][1]= F(z);
	}
	return f_value;
};
Matrix APG3(Matrix& x_0,const int &iteration,double &mu, const double &L)
{
	int n= x_0.get_row_dimension();
	Matrix f_value(iteration,2);
	Matrix x;
	x= x_0; 
	Matrix z;
	z= x;
	Matrix y;
	double A= 0.0;
	Matrix alpha(n,1) ;
	for (int i=0; i< iteration; i++){
		apg3(x,z,y,A,alpha,x_0,mu,L);
		f_value.data[i][0]= F(x);
		f_value.data[i][1]= F(z);
	}
	return f_value;
};
Matrix AdaMAPG3(Matrix &y_ini,const int &iteration,double &mu, const double &L)
{
	int k= 0;
	int n= y_ini.get_row_dimension();
	Matrix f_value(iteration,2);
	Matrix x_0;
	x_0= prox(y_ini-(1/L)*grad_f(y_ini),L);
	for (int t= 1; t< iteration; t++){
		Matrix x;
		x= x_0;
		Matrix z;
		z= x_0;
		Matrix y;
		double A= 0.0;
		Matrix alpha(n,1); 
		Matrix v;
		for (int i=0; i< Max; i++){
			apg3(x,z,y,A,alpha,x_0,mu,L);
			f_value.data[k][0]= F(x);
			f_value.data[k][1]= F(z);
			k= k+1;
			v= prox(x-(1/L)*grad_f(x),L);
			if (squarevector(v- x)<= 0.5*squarevector(x_0- y_ini)){
				x_0= x;
				y_ini= y;
				break;
			}
			else if (A>= 8*mu*mu/L){
				mu= 0.5*mu;
				break;
			}
		}
		if (k>= iteration){
			break;
		}
	}
	return f_value;
};
Matrix AdaMAPG1(Matrix &y_ini,const int &iteration,double &mu, const double &L)
{
	int num= 0;
	int k= 0;
	Matrix f_value(iteration,2);
	Matrix x_0;
	x_0= prox(y_ini-(1/L)*grad_f(y_ini),L);
	for (int t= 1; t< iteration; t++){
		Matrix x;
		x= x_0;
		Matrix y;
		y= x_0;
		Matrix y_temp;
		double theta= sqrt(mu/L);
		for (int i=0; i< Max; i++){
			y_temp= y;
			apg1(x,y,theta,mu,L);
			f_value.data[k][0]= F(x);
			f_value.data[k][1]= F(y);
			k= k+1;
			if (k>= iteration){
				break;
			}
			else if (squarevector(x- y_temp)<= 0.5*squarevector(x_0- y_ini)){
				x_0= x;
				y_ini= y_temp;
				break;
			}
			else if (8*L/mu*sqrt(2*pow(1- sqrt(mu/L),i))<= 1){
				mu= 0.5*mu;
				num= num+ 1;
				break;
			}
		}
		if (k>= iteration){
			break;
		}
	}
	cout<< "restart "<< num<< " times"<< endl;
	return f_value;
};
Matrix AdaRes(Matrix &x_0,const int &iteration,double &mu_F, const double &L)
{
	int num= 0;
	Matrix f_value(iteration,2);
	Matrix y_0;
	y_0= prox(x_0-(1/L)*grad_f(x_0),L);
	Matrix x;
	Matrix y; 
	Matrix v;
	double mu= 0.0;
	for (int t= 0; t< iteration; t++){
		double K= floor(2*sqrt(L*exp(1)/mu_F));
		double C= 16*L*squarevector(x_0- y_0)/mu_F;
		for (int k=0; k< Max; k++){
			x= y_0;
			y= x;
			double theta= 1.0;
			for (int i=0; i< K; i++){
				apg1(x,y,theta,mu,L);
				f_value.data[num][0]= F(x);
				f_value.data[num][1]= F(y);
				num= num+ 1;
				if (num>=iteration){
					break;
				}
			}
			v= prox(x-(1/L)*grad_f(x),L);
			if (num>= iteration){
				break;
			}
			else if (squarevector(v- x)> C*pow(theta/mu_F,k)){
				x_0= x;
				y_0= v;
				mu_F= 0.5*mu_F;
				break;
			}
		}
		if (num>= iteration){
			break;
		}
	}
	return f_value;
};
Matrix AdaAPG1(Matrix &x_0,const int &iteration, const double &L, bool choice)
{
	int k= 0;
	int n= x_0.get_row_dimension();
	Matrix f_value(iteration,2);
	Matrix y;
	Matrix x;
	Matrix x_temp;
	Matrix y_temp;
	double mu= 0.0;
	for (int t= 1; t< iteration; t++){
		x= x_0;
		y= x;
		double theta= 1.0;
		for (int i=0; i< Max; i++){
			x_temp= x;
			y_temp= y;
			apg1(x,y,theta,mu,L);
			f_value.data[k][0]= F(x);
			f_value.data[k][1]= F(y);
			k= k+1;
			if (k>= iteration){
				break;
			}
			else if (k>1&&choice== 0){
				if (f_value.data[k-1][0]> f_value.data[k-2][0]){
					x_0= x_temp;
					break;
				}
			}
			else if (choice== 1){
				if (innerproduct(y_temp- x,x- x_temp)> 0){
					x_0= x_temp;
					break;
				}
			}
		}
		if (k>= iteration){
			break;
		}
	}
	return f_value;
};


int main()
{
	srand (time(NULL)); //set random seed
	double L= 0.0;
	for(int i=0;i<A.row;i++)  //init the matrix
	{
		b.data[i][0]=rand() % 100/10.0;
		for(int j=0;j<A.column;j++)
		{
			A.data[i][j]=(rand() % 100/10.0);
			L= L+ A.data[i][j]*A.data[i][j];
		}
	}	
	int selection;
	int iteration;
	Matrix x_0(M,1);
	cout<<"Please enter the number of iteration:"<<endl;
	cin>>iteration;


	Matrix F_value;
	string pathnamestore;
	string p;
//	cout<<"Output?(0 or 1)"<<endl;
//	cin>>output;
	while(1==1)
	{
		cout<<"Please choose the alogrithm:"<<endl;
		cout<<"0.exit"<<endl; 
		cout<<"1.Proximal gardient"<<endl;
		cout<<"2.APG1"<<endl;
		cout<<"3.APG2"<<endl;
		cout<<"4.APG3"<<endl;
		cout<<"5.AdaMAPG3"<<endl;
		cout<<"6.AdaMAPG1"<<endl;
		cout<<"7.AdaRes"<<endl;
		cout<<"8.AdaAPG1"<<endl;
		cin>>selection;
		if(selection==0) 
		{
			break;
		}
		if(selection==1){
			F_value= PG(x_0,iteration,L);
			ofstream myfile;
			myfile.open ("PG.csv");
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<endl;
			}
			myfile.close();
		}
		if(selection==2)
		{
			double mu= 0.0;
			F_value= APG1(x_0,iteration,mu,L);
			ofstream myfile;
			myfile.open ("APG1.csv");
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==3)
		{
			double mu= 0.0;
			F_value= APG2(x_0,iteration,mu,L);
			ofstream myfile;
			myfile.open ("APG2.csv");
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}	
		if(selection==4)
		{
			double mu= 0.0;
			F_value= APG3(x_0,iteration,mu,L);
			ofstream myfile;
			myfile.open ("APG3.csv");
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==5)
		{
			double mu= 0.1*L;
			F_value= AdaMAPG3(x_0,iteration,mu,L);
			ofstream myfile;
			myfile.open ("AdaMAPG3_10per.csv");
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==6)
		{
			double est= 1e-1;
			double mu= est*L;
			F_value= AdaMAPG1(x_0,iteration,mu,L);
			ofstream myfile;
			p= to_string(int(-log10(est)));
			pathnamestore= "AdaMAPG1_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==7)
		{
			double est= 1e-1;
			double mu_F= est*L;
			F_value= AdaRes(x_0,iteration,mu_F,L);
			ofstream myfile;
			p= to_string(int(-log10(est)));
			pathnamestore= "AdaMAPG1_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==8)
		{
			F_value= AdaAPG1(x_0,iteration,L,0);
			ofstream myfile;
			myfile.open ("AdaAPG1_fun.csv");
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}	
	}
	system("pause");
	
}
