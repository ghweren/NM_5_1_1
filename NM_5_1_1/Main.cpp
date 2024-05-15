#pragma once
#include <iostream>
#include <vector>
#include "CauchySolvingGlobalError.h"
#include "Runge-Kutta_4th_32.h"
#define PI 3.14159265 

double function(double x, double y)
{

	return 2*x;
}

double eps, H, A, B, Yc, C;
using namespace std;
int main()
{
	cout << "Enter the eps\n";
	cin >> eps;
	cout << "Enter the H\n";
	cin >> H;
	cout << "Enter the A\n";
	cin >> A;
	cout << "Enter the B\n";
	cin >> B;
	cout << "Enter the C\n";
	cin >> C;
	cout << "Enter the Yc\n";
	cin >> Yc;
	try{
		CauchySolvingGlobalError solving(new Runge_Kutta_4th, function, 4);
		vector<double> result = solving.Solving(A, B, C, Yc, H, eps);
		cout << result[0] << ' ' << result[1] << ' ' << solving.exitCode() << ' ' << solving.exitEps() << ' ' << solving.exitH();
			
	}
	catch (invalid_argument e)
	{
		cout<<e.what();
	}
	system("pause");
}