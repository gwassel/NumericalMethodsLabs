#include "Header.h"
#include <cmath>
#include <iostream>
#include <fstream>

//uniform mesh - 1, Chebyshev mesh - 2
//PointVector  size = n + 1, must count in x0!
void MakeMesh(double x0, double xN, double*& PointVector, int size, int MeshType) {
	double step;
	const double PI = 3.1415;
	if (size != 0) {
		step = (xN - x0) / (size-1.0);
	}
	else {
		return;
	}

	switch (MeshType) {
	case 1:
		PointVector[0] = x0;
		for (int i = 1; i < size; i++) {
			PointVector[i] = x0 + i * step;
		}
		break;
	case 2:
		double A = (xN + x0) / 2.0;
		double B = (xN - x0) / 2.0;
		double C = 2.0 * (size + 1.0);
		std::cout << size << std::endl;
		for (int i = 0; i < size; i++) {
			PointVector[i] = A + B * cos(PI * (2.0 * i + 1.0) / C);
		}
		MergeSort(PointVector, 0, size-1, size);
	}

	double* y = new double[size];

	for (int i = 0; i < size; i++) {
		y[i] = PointVector[i] * PointVector[i];
	}

	WriteCoords("res/coord.txt",PointVector,y,size);
}


void Merge(double*& A, int first, int last, int size)
{
	int middle, start, final;
	double* mas = new double[size];
	middle = (first + last) / 2;
	start = first;
	final = middle + 1; 
	for (int j = first; j <= last; j++) 
		if ((start <= middle) && ((final > last) || (A[start] < A[final])))
		{
			mas[j] = A[start];
			start++;
		}
		else
		{
			mas[j] = A[final];
			final++;
		}
	
	for (int j = first; j <= last; j++) A[j] = mas[j];
	delete[]mas;
};

void MergeSort(double*& A, int first, int last, int size)
{
	{
		if (first < last)
		{
			MergeSort(A, first, (first + last) / 2, size); 
			MergeSort(A, (first + last) / 2 + 1, last, size); 
			Merge(A, first, last,size); 
		}
	}
};

void WriteCoords(const std::string fileNameOutput, double*& vectorX, double*& vectorY, const size_t n)
{
	std::ofstream fileOutput;
	fileOutput.open(fileNameOutput);

	for (int i = 0; i < n; ++i)
	{
		fileOutput << vectorX[i]<<" "<<vectorY[i] << "\n";
	}
}
















