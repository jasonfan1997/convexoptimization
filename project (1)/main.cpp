#include<iostream>
#include<string>
#include <cmath>
#include <time.h>
#include<vector>
#include <stdlib.h> 
//#include "MatrixXf.h"
#include <Eigen/Dense>  
#include "DENSE_DATA_READ.h"
#include <random>
#include <fstream>
#include <algorithm>
using namespace std;

int N; //define the problem
int M;
MatrixXf A;
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
	for(int i=0;i<x.get_row_dimension();i++)
	{
		result+=abs(x.data[i][0]);
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
	MatrixXf result=x;
	for(int i=0;i<x.row;i++)
	{
		if(x.data[i][0]>=0)
		result.data[i][0]=max(abs(x.data[i][0])-1/L,0.0);
		else 
		result.data[i][0]=-max(abs(x.data[i][0])-1/L,0.0);
	}
	return result;
};
MatrixXf PG(MatrixXf& x,const int &iteration,const double &L)
{
	MatrixXf f_value(iteration,1);
	for (int i=0; i< iteration; i++){
		x=prox(x-(1/L)*grad_f(x),L);
		f_value.data[i][0]= F(x);
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
MatrixXf APG1(MatrixXf& x,const int &iteration,double &theta, double &mu, const double &L)
{
	MatrixXf f_value(iteration,2);
	MatrixXf y;
	y= x;
	for (int i=0; i< iteration; i++){
		apg1(x,y,theta,mu,L);
		f_value.data[i][0]= F(x);
		f_value.data[i][1]= F(y);
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
		f_value.data[i][0]= F(x);
		f_value.data[i][1]= F(z);
	}
	return f_value;
};
MatrixXf APG3(MatrixXf& x_0,const int &iteration,double &mu, const double &L)
{
	int n= x_0.get_row_dimension();
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
		f_value.data[i][0]= F(x);
		f_value.data[i][1]= F(z);
	}
	return f_value;
};
MatrixXf AdaMAPG3(MatrixXf &y_ini,const int &iteration,double &mu, const double &L)
{
	int k= 0;
	int n= y_ini.get_row_dimension();
	MatrixXf f_value(iteration,2);
	MatrixXf x_0;
	x_0= prox(y_ini-(1/L)*grad_f(y_ini),L);
	for (int t= 1; t< iteration; t++){
		MatrixXf x;
		x= x_0;
		MatrixXf z;
		z= x_0;
		MatrixXf y;
		double A= 0.0;
		MatrixXf alpha(n,1); 
		MatrixXf v;
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
MatrixXf AdaMAPG1(MatrixXf &y_ini,const int &iteration,double &mu, const double &L)
{
	int num_1= 0;
	int num_2= 0;
	int k= 0;
	int res_1[100];
	res_1[0]= 0;
	int res_2[100];
	res_2[0]= 0;
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
			f_value.data[k][0]= F(x);
			f_value.data[k][1]= F(y);
			k= k+1;
			if (k>= iteration){
				break;
			}
			else if (squarevector(x- y_temp)<= 0.5*squarevector(x_0- y_ini)){
				x_0= x;
				y_ini= y_temp;
				num_1= num_1+ 1;
				res_1[num_1]= k;
				break;
			}
			else if (8*L/mu*sqrt(2*pow(1- sqrt(mu/L),i))<= 1){
				mu= 0.5*mu;
				num_2= num_2+ 1;
				res_2[num_2]= k;
				break;
			}
		}
		if (k>= iteration){
			break;
		}
	}
	cout<< "In all restart at"<<  endl;
	for (int i=0; i<= num_1; i++){
		cout<< res_1[i]<< "th "; 
	}
	cout<< " iteration by condition A and "<< endl;
	for (int i=0; i<= num_2; i++){
		cout<< res_2[i]<< "th "; 
	}
	cout<< " iteration by condition B"<< endl;
	return f_value;
};
MatrixXf FixRes(MatrixXf &x_0, const int &iteration, double &mu_F, const double &L) {
	int num= 0;
	MatrixXf f_value(iteration,2);
	MatrixXf x;
	x= x_0;
	MatrixXf y; 
	int K= 0;
	double mu= 0.0;
	for (int i=1; i< iteration; i++){
		K= floor(2*sqrt(L*exp(1)/mu_F));
		y= x;
		double theta= 1.0;
		for (int k=0; k< K; k++){
			apg1(x,y,theta,mu,L);
			f_value.data[num][0]= F(x);
			f_value.data[num][1]= F(y);
			num= num+ 1; 
			if (num>= iteration){
				break;
			}
		}
		if(num>= iteration){
			break;
		}
	}
	return f_value;	
};
MatrixXf Log_research(MatrixXf &x_0, const int &iteration, double &mu_F, const double &L){
	int num= 0;
	int N= ceil(4*sqrt(L*exp(1)/mu_F)*log(304*1e+5));
	int J= floor(log(N)/log(2));
	MatrixXf f_value(iteration,2);
	MatrixXf x;
	MatrixXf y; 
	double mu= 0.0;
	double K= 0.0;
	for (int t= 1; t< J; t++){
		K= pow(2,t);
		x= x_0;
		y= x;
		for (int k=0; k< Max; k++){
			double theta= 1.0;
			for (int i=0; i< K; i++){
				apg1(x,y,theta,mu,L);
				f_value.data[num][0]= F(x);
				f_value.data[num][1]= F(y);
				num= num+ 1; 
				if (num>= iteration){
					break;
				}
			}
			if (num>= iteration){
				break;
			}
			else if ((k+ 1)*K>= N){
				x_0= x;
				break;
			}
		}
		if(num>= iteration){
			break;
		}
	}
	return f_value;	
}
MatrixXf AdaRes(MatrixXf &x_0,const int &iteration,double &mu_F, const double &L, double &shr)
{
	int num= 0;
	int res_num= 0;
	int res[100]; 
	res[0]= 0;
	MatrixXf f_value(iteration,2);
	MatrixXf y_0;
	y_0= prox(x_0-(1/L)*grad_f(x_0),L);
	MatrixXf x;
	MatrixXf y; 
	MatrixXf v;
	double mu= 0.0;
	for (int t= 0; t< iteration; t++){
		double K= floor(2*sqrt(L*exp(1)/mu_F));
		double C= 16*shr*L*squarevector(x_0- y_0)/mu_F;
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
			else if (squarevector(v- x)> C*pow(theta*theta/mu_F,k)){
				x_0= x;
				y_0= v;
				mu_F= 0.5*mu_F;
				res_num= res_num+ 1;
				res[res_num]= num;
				break;
			}
			else{
				y_0= x;
			} 
		}
		if (num>= iteration){
			break;
		}
	}
	cout<< "In all restart at"<<  endl;
	for (int i=0; i<= res_num; i++){
		cout<< res[i]<< "th "; 
	}
	cout<< " iteration"<< endl;
	return f_value;
};
MatrixXf AdaAPG1(MatrixXf &x_0,const int &iteration, const double &L, int choice)
{
	int res_num= 0; 
	int k= 0;
	int n= x_0.get_row_dimension();
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
			f_value.data[k][0]= F(x);
			f_value.data[k][1]= F(y);
			k= k+1;
			if (k>= iteration){
				break;
			}
			else if (k>1&&choice== 1){
				if (f_value.data[k-1][0]> f_value.data[k-2][0]){
					x_0= x_temp;
					res_num= res_num+ 1;
					break;
				}
			}
			else if (choice== 2){
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
	cout<< "restart "<< res_num<< " times"<< endl;
	return f_value;
};
MatrixXf AdaAPG2(MatrixXf &x_0,const int &iteration, const double &L, int choice)
{
	int res_num= 0;
	int k= 0;
	int n= x_0.get_row_dimension();
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
			apg2(x,y,theta,mu,L,1);
			f_value.data[k][0]= F(x);
			f_value.data[k][1]= F(y);
			k= k+1;
			if (k>= iteration){
				break;
			}
			else if (k>1&&choice== 1){
				if (f_value.data[k-1][0]> f_value.data[k-2][0]){
					x_0= x_temp;
					res_num= res_num+ 1;
					break;
				}
			}
			else if (choice== 2){
				if (innerproduct(y_temp- x,x- x_temp)> 0){
					x_0= x_temp;
					res_num= res_num+ 1;
					break;
				}
			}
		}
		if (k>= iteration){
			break;
		}
	}
	cout<< "restart "<< res_num<< " times"<< endl;
	return f_value;
};

int main()
{
	srand (time(NULL)); //set random seed
	double L= 0.0;
	string pathname= " ";
	cout<<"Please input N: ";
	cin>>N;
	cout<<endl;
	cout<<"Please input M: ";
	cin>>M;
	cout<<endl;
	MatrixXf temp1(N,M);  
	MatrixXf temp2(N,1);
	cout<< "Please choose to generate data or read data"<< endl;
    cout<< "1. Generate;"<< endl;
    cout<< "2. Read."<< endl;
    int choice_data;
    cin>> choice_data;
    if (choice_data== 1){
    	DataGen(N, M, pathname);
    	DataRead(N,M,pathname,temp2,temp1);
	}
	else{
		cout << "Input the pathname of data you are going to use:" << endl;
    	cin >> pathname;
    	DataRead(N,M,pathname,temp2,temp1);
	}
	A=temp1;
    b=temp2;
	for(int i=0;i<A.row;i++)  //init the MatrixXf
	{
		for(int j=0;j<A.column;j++)
		{
			L= L+ A.data[i][j]*A.data[i][j];
		}
	}	
	Atranspose=A.transpose();
	AtransposeA=Atranspose*A;
	int selection;
	int iteration;
	MatrixXf x_0(M,1);
	cout<<"Please enter the number of iteration:"<<endl;
	cin>>iteration;
	MatrixXf F_value;
	string pathnamestore;
	string p;
	string s;
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
		cout<<"7.FixRes"<< endl;
		cout<<"8.Log_research"<<endl;
		cout<<"9.AdaRes"<<endl;
		cout<<"10.AdaAPG1"<<endl;
		cout<<"11.AdaAPG2"<<endl;
		cin>>selection;
		if(selection==0) 
		{
			break;
		}
		if(selection==1){
			x_0.setzero();
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
			x_0.setzero();
			double theta;
			cout<< "Pleasr choose the initial theta 1/n:"<<endl;
			int n;
			cin>> n;
			theta= 1.0/n;
			F_value= APG1(x_0,iteration,theta, mu,L);
			ofstream myfile;
			p= to_string(n);
			pathnamestore= pathname+"APG1_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==3)
		{
			double mu= 0.0;
			x_0.setzero();
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
			x_0.setzero();
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
			x_0.setzero();
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
			cout<< "Please choose the estimate coefficient of mu:"<< endl;
			double est;
			cin>> est;
			double mu= est*L;
			x_0.setzero();
			F_value= AdaMAPG1(x_0,iteration,mu,L);
			ofstream myfile;
			p= to_string(int(-log10(est)));
			pathnamestore= pathname+"AdaMAPG1_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==7)
		{
			cout<< "Please choose the estimate coefficient of mu:"<< endl;
			double est;
			cin>> est;
			double mu= est*L;
			x_0.setzero();
			F_value= FixRes(x_0,iteration,mu,L);
			ofstream myfile;
			p= to_string(int(-log10(est)));
			pathnamestore= pathname+"FixRes_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}		
		if(selection==8)
		{
			cout<< "Please choose the estimate coefficient of mu:"<< endl;
			double est;
			cin>> est;
			double mu= est*L;
			x_0.setzero();
			F_value= Log_research(x_0,iteration,mu,L);
			ofstream myfile;
			p= to_string(int(-log10(est)));
			pathnamestore= pathname+"Log_research_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				if (F_value.data[i][0]== 0){
					break;
				}
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}		
		if(selection==9)
		{
			cout<< "Please choose the estimate coefficient of mu:"<< endl;
			double est;
			cin>> est;
			double mu_F= est*L;
			cout<< "Please choose the shrink factor of C:"<< endl;
			double shr;
			cin>> shr;			
			x_0.setzero();
			F_value= AdaRes(x_0,iteration,mu_F,L,shr);
			ofstream myfile;
			p= to_string(int(-log10(est)));
			s= to_string(int(1/shr));
			pathnamestore= pathname+"AdaRes_"+p+"_"+s+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}
		if(selection==10)
		{
			cout<< "Please choose the restart scheme:"<< endl;
			cout<< "1. function scheme"<< endl;
			cout<< "2. gradient scheme"<< endl;
			int choice;
			cin>>  choice;
			x_0.setzero();
			F_value= AdaAPG1(x_0,iteration,L,choice);
			ofstream myfile;
			p= to_string(choice);
			pathnamestore= pathname+"AdaAPG1_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}	
		if(selection==11)
		{
			cout<< "Please choose the restart scheme:"<< endl;
			cout<< "1. function scheme"<< endl;
			cout<< "2. gradient scheme"<< endl;
			int choice;
			cin>>  choice;
			x_0.setzero();
			F_value= AdaAPG2(x_0,iteration,L,choice);
			ofstream myfile;
			p= to_string(choice);
			pathnamestore= pathname+"AdaAPG2_"+p+".csv";
			myfile.open (pathnamestore.c_str());
			for(int i=0;i<iteration;i++)
			{
				myfile<<i+1<<","<<F_value.data[i][0]<<","<<F_value.data[i][1]<<endl;
			}
			myfile.close();
		}			
	}
	system("pause");
	
}

