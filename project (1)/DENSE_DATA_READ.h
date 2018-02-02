#ifndef DENSE_DATA_READ_H
#define DENSE_DATA_READ_H

#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include<vector>
#include <stdlib.h> 
#include "Matrix.h"
using namespace std;

void DataGen(int m, int n, string &pathname)
{
	Matrix A(m,n);
	Matrix b(m,1);
	for(int i = 0; i < m; i++)
	{
		b.data[i][0]= rand()%10;
		for(int j = 0; j < n; j++)
		{
			A.data[i][j]=rand()%10/10.0;
		}
	}
	string m_s = to_string(m);
	string n_s = to_string(n);
	pathname = m_s + "_" + n_s +".txt";
	ofstream fout;
    fout.open(pathname.c_str());

    for(int i = 0; i < m ; i++)
    {
    	fout << b.data[i][0] << " ";
    	for(int j = 0; j < n; j++)
    	    fout << A.data[i][j] << " ";
        fout << "\n";
    }
    fout.close();
}

void DataRead(int m, int n,string pathname, Matrix &b, Matrix &M)
{
	ifstream fin;
	fin.open(pathname.c_str());
	M.row= m;
	M.column= n;
	b.row= m;
	b.column= 1;
	if (fin.is_open())
	{
	    string line;
	    int i = 0;
	    int j= 0;
	    while(getline(fin,line)){
	    	char *str = (char*)line.data();
		    char *seg;
		    j= 0;
		    seg= strtok(str, " ");
		    b.data[i][0]= atof(seg);
		    seg= strtok(NULL, " ");
			while( seg  != NULL)
		    {
				M.data[i][j]= atof(seg);
				j= j+1;
				seg= strtok(NULL, " ");
			}
			i= i+ 1;
		}	        
	}
}

#endif
