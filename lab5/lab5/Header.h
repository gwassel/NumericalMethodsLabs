#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//Functions
double f1(double x);
double f2(double x);
double f3(double x, double y);
double f4(double x, double y);
double f5(double x, double y);
double f6(double x, double y);

double f7(double x, double y);
double f8(double x, double y);

double f9(double x, double y);
double f10(double x, double y);

double f11(double x, double y);
double f12(double x, double y);

double ftest(double x, double y);

double f13(double x);
double f14(double x);
double f15(double x);

//Methods
int NewtonSystem(double*& F, double**& dF, double**& u, double*& buf, double*& x0, double*& x1, double*& x2,
    double (*f1)(double, double), double (*f2)(double, double), int size);
double BisectionMethod(double A, double B, double (*f)(double));
double NewtonOneEquation(double A, double B, double x, double (*f)(double));
double ChordMethodIteration(double A, double x, double (*f)(double));
void NewtonNonlinearSearch(double*& mesh3DX, double*& mesh3DY, int& meshSize3DX, int& meshSize3DY, double*& F,
    double**& dF, double**& u, double*& buf, double*& x0, double*& x1, double*& x2,
    double (*f1)(double, double), double (*f2)(double, double), int size);
void RelaxationMethod(double*& F, double**& dF, double*& x, int size);

//Support functions
double Derivative(double (*f)(double), double x);
double SecondDerivative(double (*f)(double), double x);
double PartialDerivativeX(double (*f)(double, double), double x, double y);
double PartialDerivativeY(double (*f)(double, double), double x, double y);
void DistributePointsToFiles(std::ofstream& f1_10, std::ofstream& f11_20, std::ofstream& f21_30, std::ofstream& f31_40, std::ofstream& fnoSols,
    double x, double y, int iters);

//Matrix functions
void InverseMatrix(double**& matrix, double**& matrixInversed, int size);
void MatrixMultVector(const double* const* matrix, const double* vector, double*& vectorResult, const size_t n);
double CubicMatrixNorm(const double* matrix, const size_t n);

//basic
void AllocateMemory(double*& mesh2D, double*& mesh3DX, double*& mesh3DY, double*& F, double**& u, double**& dF, double*& buf, double*& x0,
    double*& x1, double*& x2, int Fsize, int& meshSize2D, int& meshSize3DX, int& meshSize3DY);
void CalculateOneEquation(double*& table, int& size, double A, double B, double (*f)(double));

//Support functions
int HasSolution2D(double xi, double xj, double (*f)(double));
int HasSolution3D(double xi, double yi, double xj, double yj, double (*f)(double, double));
void MakeMesh2D(double A, double B, double*& table, int& size);
void PrintTable(double*& table, int& size);
void CheckForEqualElements(std::vector<double>& sols);
void CheckForEqualElements3D(std::vector<double>& solsNX, std::vector<double>& solsNY, std::vector<int>& iters);