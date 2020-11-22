#include "Header.hpp"


void ReadInit()
{

}

void AllocateMemory()
{
    
}

void ReadData()
{
    
}

void WriteData()
{
    
}

void FreeMemory()
{

}


void WriteCoords(const std::string fileNameOutput, Grid &grid)
{
	std::ofstream fileOutput;
	fileOutput.open(fileNameOutput);

	for (int i = 0; i < grid.length; ++i)
	{
		fileOutput << grid.points[i].x << " " << grid.points[i].y << "\n";
	}
}


Point::Point(double x_ /*= 1*/, double y_ /*= 1*/)
{
    x = x_;
    y = y_;
}


Grid::Grid(size_t length_ /*= 0*/)
{
    length = length_;
    points = new Point[length];
}
Grid::~Grid()
{
    delete[] points;
}
Grid::Grid(Grid&& other)
{
    points = other.points;
    other.points = nullptr;
    length = other.length;
    other.length = 0;
}


Polynomial::Polynomial(size_t length_ /*= 0*/)
{
    length = length_;
    coefficents = new double[length];
    for(int i = 1; i < length; ++i)
    {
        coefficents[i] = 0;
    }
    if(length)
    {
        coefficents[0] = 1;
    }
}
Polynomial::~Polynomial()
{
    delete[] coefficents;
}
Polynomial::Polynomial(Polynomial&& other)
{
    coefficents = other.coefficents;
    other.coefficents = nullptr;
    length = other.length;
    other.length = 0;
}
double Polynomial::eval(double x)
{
    double sum = coefficents[0];

    for(int i = 1; i < length; ++i)
    {
        sum += coefficents[i] * pow(x, i);
    }

    return sum;
}

Basis::Basis(size_t length_ /*= 0*/)
{
    length = length_;
    polynomials = new Polynomial[length];

    for(int i = 0; i < length; ++i)
    {
        polynomials[i] = Polynomial(length);
    }        

}
Basis::~Basis()
{
    delete[] polynomials;
}
Basis::Basis(Basis&& other)
{
    polynomials = other.polynomials;
    other.polynomials = nullptr;
    length = other.length;
    other.length = 0;
}

TridiagonalMatrix::TridiagonalMatrix(size_t n_)
{
    n = n_;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    d = new double[n];

    x = new double[n];

    alpha = new double[n];
    beta = new double[n];
}
TridiagonalMatrix::TridiagonalMatrix(size_t n_, double* a_, double* b_, double* c_, double* d_)
{
    n = n_;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    d = new double[n];

    x = new double[n];

    alpha = new double[n];
    beta = new double[n];

    for(int i = 0; i < n; ++i)
    {
        a[i] = a_[i];
        b[i] = b_[i];
        c[i] = c_[i];
        d[i] = d_[i];
    }
}
TridiagonalMatrix::~TridiagonalMatrix()
{
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    delete[] x;

    delete[] alpha;
    delete[] beta;
}
TridiagonalMatrix::TridiagonalMatrix(TridiagonalMatrix&& other)
{
    n = other.n;
    a = other.a;
    b = other.b;
    c = other.c;
    d = other.d;

    x = other.x;

    alpha = other.alpha;
    beta = other.beta;

    other.n = 0;
    other.a = nullptr;
    other.b = nullptr;
    other.c = nullptr;
    other.d = nullptr;
    other.alpha = nullptr;
    other.beta = nullptr;
    other.x = nullptr;
}
void TridiagonalMatrix::run()
{
    double denominator = 0;

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];

    for(int i = 1; i < n - 1; ++i)
    {
        denominator = 1 / (b[i] - a[i] * alpha[i - 1]);
        alpha[i] = c[i] * denominator;
        beta[i] = (d[i] + a[i] * beta[i - 1]) * denominator;
    }

    x[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 2]) / (b[n - 1] - a[n - 1] * alpha[n - 2]);

    for(int i = n - 2; i >= 0; --i)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}

void test(Polynomial &p1, Grid &grid, std::string label)
{
    std::cout << label << "\n";
    double err1 = CheckOnGridPoints(p1, grid);
    std::cout << "err:" << err1 << "\n\n";
}
