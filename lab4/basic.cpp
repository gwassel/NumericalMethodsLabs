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

void test(Polynomial &p1, Grid &grid, std::string label)
{
    std::cout << label << "\n";
    double err1 = CheckOnGridPoints(p1, grid);
    std::cout << "err:" << err1 << "\n\n";
}