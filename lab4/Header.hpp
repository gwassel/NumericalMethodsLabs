#include <iostream>


struct Point
{
    double x = 0;
    double y = 0;

    Point(double x_=0, double y_=0)
    {
        x = x_;
        y = y_;
    }
};

struct Grid
{
    Point* points = nullptr;
    size_t length;

    Grid(size_t length_=0)
    {
        length = length_;
        points = new Point[length];
    }

    ~Grid()
    {
        delete[] points;
    }

    Grid(const Grid& other) = delete;
    Grid(Grid&& other)
    {
        points = other.points;
        other.points = nullptr;
    }
    Grid& operator=(const Grid & other) = delete;
    Grid& operator=(Grid&& other) = delete;
};

struct Polynomial
{
    double* coefficents;
    size_t length = 0;

    Polynomial(size_t degree_=0)
    {
        length = degree_ + 1;
        coefficents = new double[length];
    }

    ~Polynomial()
    {
        delete[] coefficents;
    }

    Polynomial(const Polynomial & other) = delete;
    Polynomial(Polynomial&& other)
    {
        coefficents = other.coefficents;
        other.coefficents = nullptr;
    }

    Polynomial& operator=(const Polynomial & other) = delete;
    Polynomial& operator=(Polynomial&& other) = delete;
};

struct Basis
{
    Polynomial* polynomials = nullptr;
    size_t length = 0;

    Basis(size_t length_=0)
    {
        length = length_;
        polynomials = new Polynomial[length];
    }

    ~Basis()
    {
        delete[] polynomials;
    }

    Basis(const Basis & other) = delete;
    
    Basis(Basis&& other)
    {
        polynomials = other.polynomials;
        other.polynomials = nullptr;
    }

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


