
#include<iostream>
#include<string>
#include <cmath>
#include <time.h>
#include<vector>
#include <stdlib.h> 
//#include "Matrix.h"
#include <random>
#include <Eigen/Dense>  
#include <fstream>
#include <algorithm>
using namespace Eigen;  
using namespace std;

int N; //define the problem
int M;
MatrixXf A;
MatrixXf Atranspose;
MatrixXf AtransposeA;
MatrixXf b;
int Max= 100000000;
double squarevector(const MatrixXf &M) //taking square of the MatrixXf
{
	/*
	double result=0;
	for(int i=0;i<M.rows();i++)
	{
		result+=M(i,0)*M(i,0);
	}
	*/
	return M.squaredNorm();
};
double innerproduct(const MatrixXf &M,const MatrixXf &N) //taking inner product of the MatrixXf
{
	/*
	double result=0;
	for(int i=0;i<M.rows();i++)
	{
		result+=M(i,0)*N(i,0);
	}
	return result;
	*/
	return (M.transpose()*N).sum();
};

double f(MatrixXf &x)
{
	return 0.5*squarevector(A*x-b);
};


double g(const MatrixXf &x)
{
	double result=0;
	for(int i=0;i<x.rows();i++)
	{
		result+=abs(x(i,0));
	}
	return result;
};
double F(MatrixXf &x)
{
	return f(x)+g(x);
};

MatrixXf grad_f(MatrixXf &x)
{
	return AtransposeA*x-Atranspose*b;
};
MatrixXf prox(const MatrixXf& x,const double &L)
{
	MatrixXf result;
	result=x;
	for(int i=0;i<x.rows();i++)
	{
		if(x(i,0)>=0)
		result(i,0)=max(abs(x(i,0))-1/L,0.0);
		else 
		result(i,0)=-max(abs(x(i,0))-1/L,0.0);
	}
	return result;
};
MatrixXf PG(MatrixXf& x,const int &iteration,const double &L)
{
	MatrixXf f_value(iteration,1);
	for (int i=0; i< iteration; i++){ 
		x=prox(x-(1/L)*grad_f(x),L);	
		f_value(i,0)= f(x);
	}
	return f_value;	
};
void apg1(MatrixXf &x,MatrixXf &y,double &theta,double &mu,const double &L)
{
	MatrixXf temp;
	temp =x;
	x=prox(y-(1/L)*grad_f(y),L);
	double temptheta=theta;
	theta=(-pow(theta,2)+ mu/L+sqrt(pow(pow(theta,2)- mu/L,2)+4*pow(theta,2)))/2;
	double beta=temptheta*(1-temptheta)/(temptheta*temptheta+theta);
	y=x+beta*(x-temp);
}
void apg2(MatrixXf &x,MatrixXf &z,double &theta,double &mu, const double &L, bool choice)
{
	MatrixXf y;
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
void apg3(MatrixXf &x,MatrixXf &z,MatrixXf &y,double &A,MatrixXf &alpha,MatrixXf &x0,double &mu,const double &L)
{
	double a=(2*(1+ mu*A)/L+sqrt(4*pow((1+ mu*A)/L,2)+8*(1+ mu*A)*A/L))/2;
	y=(A/(A+a))*x+(a/(A+a))*z;
	x=prox(y-(1/L)*grad_f(y),L);
	alpha=alpha+a*grad_f(x);
	for(int i=0;i<x.rows();i++)
	{
		if((x0(i,0)-alpha(i,0))>A)
		{
			z(i,0)=x0(i,0)-alpha(i,0)-A;
		}
		else if(abs(x0(i,0)-alpha(i,0))<A)
		{
			z(i,0)=0;
		}
		else if((x0(i,0)-alpha(i,0))<-A)
		{
			z(i,0)=x0(i,0)-alpha(i,0)+A;
		}
	}
	A= A+ a;
}
MatrixXf APG1(MatrixXf& x,const int &iteration,double &mu, const double &L)
{
	MatrixXf f_value(iteration,2);
	MatrixXf y;
	y= x;
	double theta= 1.0;
	for (int i=0; i< iteration; i++){
		apg1(x,y,theta,mu,L);
		f_value(i,0)= F(x);
		f_value(i,1)= F(y);
	}
	return f_value;
};
MatrixXf APG2(MatrixXf& x,const int &iteration,double &mu, const double &L)
{
	MatrixXf f_value(iteration,2);
	MatrixXf z;
	z= x;
	double theta= 1.0;
	for (int i=0; i< iteration; i++){
		apg2(x,z,theta,mu,L,0);
		f_value(i,0)= F(x);
		f_value(i,1)= F(z);
	}
	return f_value;
};
MatrixXf APG3(MatrixXf& x_0,const int &iteration,double &mu, const double &L)
{
	int n= x_0.rows();
	MatrixXf f_value(iteration,2);
	MatrixXf x;
	x= x_0; 
	MatrixXf z;
	z= x;
	MatrixXf y;
	double A= 0.0;
	MatrixXf alpha(n,1) ;
	for (int i=0; i< iteration; i++){
		apg3(x,z,y,A,alpha,x_0,mu,L);
		f_value(i,0)= F(x);
		f_value(i,1)= F(z);
	}
	return f_value;
};
MatrixXf AdaMAPG3(MatrixXf &y_ini,const int &iteration,double &mu, const double &L)
{
	int k= 0;
	int n= y_ini.rows();
	MatrixXf f_value(iteration,2);
	MatrixXf x_0;
	x_0= prox(y_ini-(1/L)*grad_f(y_ini),L);
	
	for (int t= 1; t< iteration; t++){
		MatrixXf x;
		x= x_0;
		MatrixXf z;
		z= x_0;
		MatrixXf y=y_ini;
		double A= 0.0;
		MatrixXf alpha=MatrixXf::Zero(n,1); 
		MatrixXf v;
		for (int i=0; i< Max; i++){
			
			apg3(x,z,y,A,alpha,x_0,mu,L);
			f_value(k,0)= F(x);
			f_value(k,1)= F(z);
			k= k+1;
			cout<<k<<endl;
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
		cout<<"temp2"<<endl;
	}
	return f_value;
};
MatrixXf AdaMAPG1(MatrixXf &y_ini,const int &iteration,double &mu, const double &L)
{
	int num= 0;
	int k= 0;
	MatrixXf f_value(iteration,2);
	MatrixXf x_0;
	x_0= prox(y_ini-(1/L)*grad_f(y_ini),L);
	for (int t= 1; t< iteration; t++){
		MatrixXf x;
		x= x_0;
		MatrixXf y;
		y= x_0;
		MatrixXf y_temp;
		double theta= sqrt(mu/L);
		for (int i=0; i< Max; i++){
			y_temp= y;
			apg1(x,y,theta,mu,L);
			f_value(k,0)= F(x);
			f_value(k,1)= F(y);
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
MatrixXf AdaRes(MatrixXf &x_0,const int &iteration,double &mu_F, const double &L)
{
	int num= 0;
	MatrixXf f_value(iteration,2);
	MatrixXf y_0;
	y_0= prox(x_0-(1/L)*grad_f(x_0),L);
	MatrixXf x;
	MatrixXf y; 
	MatrixXf v;
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
				f_value(num,0)= F(x);
				f_value(num,1)= F(y);
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
MatrixXf AdaAPG1(MatrixXf &x_0,const int &iteration, const double &L, bool choice)
{
	int k= 0;
	int n= x_0.rows();
	MatrixXf f_value(iteration,2);
	MatrixXf y;
	MatrixXf x;
	MatrixXf x_temp;
	MatrixXf y_temp;
	double mu= 0.0;
	for (int t= 1; t< iteration; t++){
		x= x_0;
		y= x;
		double theta= 1.0;
		for (int i=0; i< Max; i++){
			x_temp= x;
			y_temp= y;
			apg1(x,y,theta,mu,L);
			f_value(k,0)= F(x);
			f_value(k,1)= F(y);
			k= k+1;
			if (k>= iteration){
				break;
			}
			else if (k>1&&choice== 0){
				if (f_value(k-1,0)> f_value(k-2,0)){
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
	double L= 1000;
	cout<<"Please input N: ";
	cin>>N;
	cout<<endl;
	cout<<"Please input M: ";
	cin>>M;
	cout<<endl;
	MatrixXf temp1(N,M);  
	MatrixXf temp2(N,1);
	//MatrixXf temp3(M,M);
	A=temp1;
    b=temp2;
    
	cout<<"Please select the problem(0 for random generated A and b,1 for user input A and b,2 for sparse"<<endl;
	int problem;
	cin>>problem;
	if(problem==0)
	{	
		for(int i=0;i<A.rows();i++)  //init the MatrixXf
		{
			
		b(i,0)=rand() % 100/10.0;
			for(int j=0;j<A.cols();j++)
			{
				if(i==j)
					A(i,j)=rand() % 100/10.0;
			}
		}
	}
	else if (problem==1)
	{
		cout<<"Please input MatrixXf A and b:"<<endl;
		for(int i=0;i<N;i++)  //init the MatrixXf
		{
			cin>>b(i,0);
			for(int j=0;j<M;j++)
			{
					cin>>A(i,j);
			}
		}
	}
	else
	{
		cout<<"Please input MatrixXf A and b in sprase form:(enter end at the end)"<<endl;
		string temp;
		cin>>temp;
		int i=0;
		while(i<N)
		{
			
			b(i,0)=stod(temp);
			cin>>temp;
			while(temp.find(":")!= std::string::npos && temp!="end")
			{
				int j;
				for(j =0;j<temp.size();j++)
				{
					if(temp[j]==':')
						break; 
				}
				A(i,stoi(temp.substr(0,j)))=stod(temp.substr(j+1,temp.size()-j));
				cin>>temp;
			}
			i++;
		}
	}
//	cout<<"temp"<<endl;
	Atranspose=A.transpose();
	AtransposeA=Atranspose*A;
	int selection;
	int iteration;
	cout<<"Please enter the number of iteration:"<<endl;
	cin>>iteration;


	
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
		MatrixXf F_value;
		cin>>selection;
		MatrixXf x_0=MatrixXf::Zero(M,1);
		//cout<<x_0<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
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
				myfile<<i+1<<","<<F_value(i,0)<<","<<F_value(i,1)<<endl;
			}
			myfile.close();
		}	
	}
	system("pause");
	
}

