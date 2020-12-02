#include "Header.hpp"


void ReadInit()
{

}

void ReadData()
{
    
}

void WriteData()
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
    fileOutput.close();
}

void WritePolynomial(const std::string fileNameOutput, Polynomial &p1, Grid &grid)
{
    std::ofstream fileOutput;
	fileOutput.open(fileNameOutput);

    fileOutput << "{";
	for (int i = 0; i < grid.length; ++i)
	{
		fileOutput << std::fixed << "{" << grid.points[i].x << ", " << p1.eval(grid.points[i].x) << "},\n";
	}
    fileOutput << "}";
    fileOutput.close();
}

void WriteSpline(const std::string fileNameOutput, Spline &spline)
{
    std::ofstream fileOutput;
	fileOutput.open(fileNameOutput);
    size_t stepSize = 25;
    double x0 = spline.x[0];
    double x1 = spline.x[1];
    double step = (x1 - x0) / stepSize;
    double xcurrent;
    for(int i = 0; i < spline.n; ++i)
    {   
        double x0 = spline.x[i];
        double x1 = spline.x[i+1];
        double step = (x1 - x0) / stepSize;
        xcurrent = x0;
        for(int j = 0; j <= stepSize; ++j)
        {
            fileOutput << xcurrent << " " << spline.a[i] + spline.b[i] * (xcurrent - x0) + spline.c[i] * pow((xcurrent - x0), 2) + spline.d[i] * pow(xcurrent - x0,3) << "\n";
            xcurrent += step;
        }
    }
    xcurrent = spline.x[spline.n];
    x0 = spline.x[spline.n-1];
    fileOutput << xcurrent << " " << spline.a[spline.n-1] + spline.b[spline.n-1] * (xcurrent - x0) + spline.c[spline.n-1] * pow((xcurrent - x0), 2) + spline.d[spline.n-1] * pow(xcurrent - x0,3) << "\n";
    fileOutput.close();
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


TridiagonalMatrix::TridiagonalMatrix(double* &h, double* &g, size_t n_)// n - size of h
{
    n = n_ - 1;
    
    a = new double[n];
    b = new double[n];
    c = new double[n];
    d = new double[n];

    x = new double[n];

    alpha = new double[n];
    beta = new double[n];


    for(int i = 1; i < n; ++i)
    {
        a[i] = h[i];
    }

    for(int i = 0; i < n - 1; ++i)
    {
        c[i] = h[i + 1];
    }

    for(int i = 0; i < n; ++i)
    {
        b[i] = 2 * (h[i] + h[i + 1]);
    }

    for(int i = 0; i < n; ++i)
    {
        d[i] = 3 * (g[i + 1] - g[i]);
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
void TridiagonalMatrix::eval()
{
    // for(int i = 0; i < n; ++i)
    // {
    //     b[i] = -b[i];
    //     d[i] = -d[i];
    // }
}
void TridiagonalMatrix::run()
{
    double denominator = 0;

    for(int i = 1; i < n - 1; ++i)
    {
        denominator = a[i]/b[i-1];
        b[i] = b[i] - denominator * c[i-1];
        d[i] = d[i] - denominator * d[i-1];
    }

    x[n-2] = d[n-2]/b[n-2];

    for(int i = n - 3; i >= 0; --i)
    {
        x[i] = (d[i] - c[i]*x[i+1])/b[i];
    }
}


Spline::Spline()
{
    n = 0;

    double* a = nullptr;
    double* b = nullptr;
    double* c = nullptr;
    double* d = nullptr;

    double* h = nullptr;
    double* g = nullptr;

    double* x = nullptr;
}
Spline::Spline(Grid& grid)
{
    n = grid.length - 1;

    a = new double[n];
    b = new double[n];
    c = new double[n];
    d = new double[n];

    h = new double[n];
    g = new double[n];

    x = new double[n+1];


    for(int i = 0; i < n; ++i)
    {
        x[i] = grid.points[i].x;
        h[i] = grid.points[i + 1].x - grid.points[i].x;
        g[i] = (grid.points[i + 1].y - grid.points[i].y) / h[i];
        a[i] = grid.points[i].y;
    }
    x[n]  = grid.points[n].x;
}

Spline::~Spline()
{
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    delete[] h;
    delete[] g;

    delete[] x;
}

Spline::Spline(Spline&& other)
{
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
}
void Spline::RecountCoefficents(TridiagonalMatrix& matrix)
{
    b[0] = g[0] - matrix.x[0] * h[0] / 3;
    d[0] = matrix.x[0] / 3 / h[0];
    c[0] = 0;

    for(int i = 1; i < n; ++i)
    {
        c[i] = matrix.x[i - 1];
    }

    for(int i = 1; i < n - 1; ++i)
    {
        b[i] = g[i] - (c[i+1] + 2.0 * c[i]) / 3.0 * h[i];
        d[i] = (c[i+1] - c[i]) / 3.0 / h[i];
    }

    b[n - 1] = g[n - 1] - 2 * matrix.x[n - 2] / 3 * h[n - 1];

    d[n - 1] = -1.0 / 3.0 / h[n - 1] * c[n - 1];
}
double Spline::eval(double x_)
{
    int i = 0;
    while(x[i] < x_)
    {
        x_ = x[i];
        ++i;
    }
    --i;
    return a[i] + b[i] * (x_ - x[i]) + c[i] * pow((x_ - x[i]), 2) + d[i] * pow((x_ - x[i]), 3);
}

void test(Polynomial &p1, Grid &grid, std::string label)
{
    std::cout << label << "\n";
    double err1 = CheckOnGridPoints(p1, grid);
    std::cout << "err:" << err1 << "\n\n";
}

double CountError(Polynomial &p, Grid &testGrid, double (*f)(double), double leftBorder, double rightBorder)
{
    double error = 0;
    double errIterations = 0;
    
    MakeMesh(leftBorder, rightBorder, testGrid, 1, f);

    for(int i = 0; i < testGrid.length; ++i)
    {
        errIterations = fabs(p.eval(testGrid.points[i].x) - testGrid.points[i].y);
        if( errIterations > error )
        {
            error = errIterations;
        }
    }    

    return error;
}
double CountError(Spline &spline, double (*f)(double))
{
    size_t stepSize = 25;
    double x0 = spline.x[0];
    double x1 = spline.x[1];
    double step = (x1 - x0) / stepSize;
    double xcurrent;
    double err = 0;
    double errIter = 0;
    for(int i = 0; i < spline.n; ++i)
    {   
        double x0 = spline.x[i];
        double x1 = spline.x[i+1];
        double step = (x1 - x0) / stepSize;
        xcurrent = x0;
        for(int j = 0; j <= stepSize; ++j)
        {
            errIter = fabs(spline.a[i] + spline.b[i] * (xcurrent - x0) + spline.c[i] * pow((xcurrent - x0), 2) + spline.d[i] * pow(xcurrent - x0,3) - f(xcurrent));
            if(errIter > err)
            {
                err = errIter;
            }
            xcurrent += step;
        }
    }
    return err;
}

