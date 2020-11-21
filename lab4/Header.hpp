#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>

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
    Polynomial& operator=(Polynomial&& other) = delete;
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

void ReadInit();
void AllocateMemory();
void ReadData();

void MakeUniformMesh();
void MakeChebishevMesh();//name may be incorrect

void LagrangeInterpolate(Grid& grid, Basis &basis, Polynomial &pLagrange, Polynomial &pBuffer);
void SplineInterpolate();

void WriteData();
void FreeMemory();

void MakeMesh(double x0, double xN, Grid &grid, int MeshType);
void MergeSort(double*& A, int first, int last, int size);
void Merge(double*& A, int first, int last, int size);
void WriteCoords(const std::string fileNameOutput, Grid &grid);
