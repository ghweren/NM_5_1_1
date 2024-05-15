#pragma once
#include "SolvingMethod.h"
#include <vector>
#include <cmath>
using namespace std;

class CauchySolvingGlobalError
{
private:
	
	double(*function)(double,double);
	
	SolvingMethod *solvingMethod;
	
	size_t lcod;

	double macheps, H, eps, maxEps, hMin, A,B,s;

	double eps_order(double eps);

	void machepsCalculation()
	{
		double R = 1;
		while (1 + R > 1)
			R /= 2;
		macheps = R * 2;
	}

	void hMimCalculation()
	{
		hMin = macheps * std::fmax(abs(A), std::fmax(abs(B), DBL_MIN));
	}

	double HCalculation(double H, double y1, double y2)
	{
		return H / 2 * pow((pow(2, s) - 1)*eps / abs(y1 - y2), 1 / s);
	}
	

public:
	CauchySolvingGlobalError(SolvingMethod *solvingMethod, double(*function)(double, double), size_t order)
	{
		this->solvingMethod = solvingMethod;
		this->function = function;
		machepsCalculation();
		this->s = order;
	}

	
	vector<double> Solving(double A,double B, double C, double Yc,double H,double eps);

	size_t exitCode()
	{
		return lcod;
	}

	double exitH()
	{
		return H;
	}

	double exitEps()
	{
		return eps;
	}


};
