#include "CauchySolvingGlobalError.h"
#include <iostream>

vector<double> CauchySolvingGlobalError::Solving(double A, double B, double C, double Yc, double H, double eps)
{
	vector<double> result(2);
	double prev_error;
	double curr_error=0;
	maxEps = eps;
	hMimCalculation();
	H = (B - A) / ceil((B - A) / H);
	this->H = H;
	if (H <= hMin)
		throw new invalid_argument("H<=hMin");
	do
	{
		prev_error = curr_error;
		if (C == A)
		{
			double prev_point = A;
			double prev_result = Yc,prev_estimationResult=Yc;
			double int_result;
			double estimationResult;
			while (abs(prev_point-B)>eps)
			{
				int_result = solvingMethod->Calculation(prev_result, prev_point, H, function);
				estimationResult = solvingMethod->Calculation(prev_estimationResult, prev_point, H / 2, function);
				estimationResult = solvingMethod->Calculation(estimationResult, prev_point + H / 2, H / 2, function);
				prev_point += H;
				prev_result = int_result;
				prev_estimationResult = estimationResult;
			}
			curr_error = abs(int_result - estimationResult);
			this->H = H;
			H = HCalculation(H,int_result,estimationResult);
			result[0] = B;
			result[1] = int_result;
		}
		else
		{
			double prev_point = B;
			double prev_result = Yc, prev_estimationResult = Yc;
			double int_result;
			double estimationResult;
			while (abs(prev_point-A)>eps)
			{
				int_result = solvingMethod->Calculation(prev_result, prev_point, -H, function);
				estimationResult = solvingMethod->Calculation(prev_estimationResult, prev_point, -H / 2, function);
				estimationResult = solvingMethod->Calculation(estimationResult, prev_point - H / 2, -H / 2, function);
				prev_point -= H;
				prev_result = int_result;
				prev_estimationResult = estimationResult;
			}
			curr_error = abs(int_result - estimationResult);
			this->H = H;
			H = HCalculation(H, int_result, estimationResult);
			result[0] = A;
			result[1] = int_result;
		}
	} while (H > hMin && curr_error > maxEps && abs(curr_error - prev_error) > maxEps);
	eps = curr_error;
	if (abs(curr_error - prev_error) <= maxEps)
		lcod = 1;
	else
	if (curr_error <= maxEps)
		lcod = 0;
	else
		lcod = 2;
	return result;
}
double CauchySolvingGlobalError::eps_order(double eps)
{
	if (eps < 1)
	{
		int order = 0;
		while (eps < 1)
		{
			eps *= 10;
			order--;
		}
		return order;
	}
	else
	{
		int order = 0;
		while (eps>1)
		{
			eps /= 10;
			order++;
		}
		return order;
	}
}