#pragma once
class SolvingMethod abstract
{
public:
	virtual double Calculation(double Y0, double X0, double h, double(*function)(double, double)) const = 0;
};