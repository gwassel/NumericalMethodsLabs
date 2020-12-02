#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "stdio.h"
#include "string.h"
#include "math.h"

struct Point
{
    double x = 0;
    double y = 0;

    Point(double x_=0, double y_=0);
};

struct Grid
{
    Point* points = nullptr;
    size_t length = 0;

    Grid(size_t length_=0);
    ~Grid();
    Grid(const Grid& other) = delete;
    Grid(Grid&& other);
    Grid& operator=(const Grid & other) = delete;
    Grid& operator=(Grid&& other) = delete;
};

struct Polynomial
{
    double* coefficents = nullptr;
    size_t length = 0;

    Polynomial(size_t length_=0);
    ~Polynomial();
    Polynomial(const Polynomial & other) = delete;
    Polynomial(Polynomial&& other);
    Polynomial& operator=(const Polynomial & other) = delete;
    Polynomial& operator=(Polynomial&& other)
    {
        if (this == &other)
            return *this;

        delete[] coefficents;
        coefficents = other.coefficents;
        length = other.length;
        other.coefficents = nullptr;
        other.length = 0;
        return *this;
    }

    double eval(double x);
};

struct Basis
{
    Polynomial* polynomials = nullptr;
    size_t length = 0;

    Basis(size_t length_=0);
    ~Basis();
    Basis(const Basis & other) = delete;
    Basis(Basis&& other);
    Basis& operator=(const Basis& other) = delete;
    Basis& operator=(Basis&& other) = delete;
};

struct TridiagonalMatrix
{
    size_t n = 0;

    double* a = nullptr; //subdiagonal
    double* b = nullptr; //diagonal
    double* c = nullptr; //superdiagonal
    double* d = nullptr; //right side vector

    double* x = nullptr;

    double* alpha = nullptr;
    double* beta = nullptr;

    TridiagonalMatrix(double* &h, double* &g, size_t n);
    ~TridiagonalMatrix();
    TridiagonalMatrix(const TridiagonalMatrix& other) = delete;
    TridiagonalMatrix(TridiagonalMatrix&& other);
    TridiagonalMatrix& operator=(const TridiagonalMatrix& other) = delete;
    TridiagonalMatrix& operator=(TridiagonalMatrix&& other) = delete;

    void eval();
    void run();
};

struct Spline
{
    size_t n = 0;

    double x0 = 0;

    double* a = nullptr;
    double* b = nullptr;
    double* c = nullptr;
    double* d = nullptr;

    double* x = nullptr;

    double* h = nullptr;
    double* g = nullptr;
    
    void RecountCoefficents(TridiagonalMatrix& matrix);

    Spline();
    Spline(Grid &grid);
    ~Spline();
    Spline(const Spline& other) = delete;
    Spline(Spline&& other);
    Spline& operator=(const Spline& other) = delete;
    Spline& operator=(Spline&& other)
    {
        if (this == &other)
            return *this;

        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
        delete[] h;
        delete[] g;
        delete[] x;

        n = other.n;
        a = other.a;
        b = other.b;
        c = other.c;
        d = other.d;
        h = other.h;
        g = other.g;
        x = other.x;

        other.n = 0;
        other.a = nullptr;
        other.b = nullptr;
        other.c = nullptr;
        other.d = nullptr;
        other.h = nullptr;
        other.g = nullptr;
        other.x = nullptr;
        
        return *this;
    }
    double eval(double x);
};

void ReadInit();
void AllocateMemory();
void ReadData();

void LagrangeInterpolate(Grid &grid, Basis &basis, Polynomial &pLagrange, Polynomial &pBuffer, Polynomial &monomial, std::string label = "");
void SplineInterpolate();
double CheckOnGridPoints(Polynomial &p1, Grid &grid);

void WriteData();
void FreeMemory();

void MakeMesh(double x0, double xN, Grid &grid, int MeshType, double (*f)(double));
void MergeSort(double*& A, int first, int last, int size);
void Merge(double*& A, int first, int last, int size);
void WriteCoords(const std::string fileNameOutput, Grid &grid);
void WriteSpline(const std::string fileNameOutput, Spline &spline);
void WritePolynomial(const std::string fileNameOutput, Polynomial &p1, Grid &grid);

void test(Polynomial &p1, Grid &grid, std::string label);

//Count functions
double f0(double);
double f1(double);
double f2(double);
double f3(double);
double f4(double);
double f5(double);
double f6(double);
double f7(double);
double f8(double);
double f9(double);

int MakeSpline(Grid &grid);



void LagrangeInterpolate(Grid &initGrid, Grid &lagrange);
void LagrangeOut(Grid &lagrangeGrid, std::string fileName);
double CountError(Grid &grid, double (*f)(double));
double CountError(Spline &spline, double (*f)(double));
