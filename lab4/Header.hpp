#include <iostream>


struct Point
{
    double x = 0;
    double y = 0;

    Point(double x_, double y_)
    {
        x = x_;
        y = y_;
    }
};

struct Grid
{
    Point* points = nullptr;
    size_t length;

    Grid(size_t length_)
    {
        length = length_;
        points = new Point[length];
    }

    ~Grid()
    {
        delete[] points;
    }
};

struct Polynomial
{
    double* coefficents = nullptr;
    size_t length = 0;

    Polynomial(size_t length_)
    {
        length = length_;
        coefficents = new double[length];
    }

    ~Polynomial()
    {
        delete[] coefficents;
    }
};

void ReadInit();
void AllocateMemory();
void ReadData();

void MakeUniformMesh();
void MakeChebishevMesh();//name may be incorrect

void LagrangeInterpolate();
void SplineInterpolate();

void WriteData();
void FreeMemory();


