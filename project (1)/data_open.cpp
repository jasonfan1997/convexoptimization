#include "DENSE_DATA_READ.h"

using namespace std;

int main()
{
//------------------Construction of Linear System------------------
	string pathname;
    cout << "Input the pathname of data you are going to use:" << endl;
    cin >> pathname;
    int m;
    cout << "Input the row number of your data:" << endl;
    cin >> m;

    int n;
    cout << "Input the column number of your data:" << endl;
    cin >> n;
/*	std::vector<string> Data1, Data2;
	Data1= DataGen<double, int>(1000, 10000, 5000, pathname, 0);
	Data2= DataGen<double, int>(100, 200, 100, pathname, 0);
	Parameters<double,int> data;
	data= DataRead<double, int>(2,8,pathname, 1);
	MATRIX<double, int> A;
	A= data.M; 
	vector<double> c2= data.c;
	for (int i=0; i< c2.size(); i++){
		cout<< c2[i]<< "/s";
	}
	vector<double> d2= data.d;
	for (int i=0; i< d2.size(); i++){
		cout<< d2[i]<< "/s";
	}*/
	DataGen(m, n, pathname);
	Matrix b(m,1);
	Matrix M(m,n);
	DataRead(m,n,pathname,b,M);
	cout<< "read done!"<< endl;
	for (int i=0; i< m; i++){
		cout<< "the "<< i<< "th number of c is "<< b.data[i][0]<< endl; 
	}
	for (int i=0; i< m; i++){
		for (int j=0; j< n; j++){
			cout<< "the "<< j << "th number of"<< i<< "th row of M is"<< M.data[i][j]<< endl;
			system("pause");
		}
	}
 //	LS<double, int> lp(2, 8, pathname);
	return 0;
}
