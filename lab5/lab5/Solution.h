#pragma once
#include "Header.h"
class Solution
{
private:
	std::vector<double> X;
	std::vector<double> Y;
	void Add(double x, double y) {
		X.push_back(x);
		Y.push_back(y);
	};
	bool IsEqual(double x, double y) {
		double eps = 1e-8;};
public:
	double getXi(int i) { return X[i]; };
	double getYi(int i) { return Y[i]; };

};

