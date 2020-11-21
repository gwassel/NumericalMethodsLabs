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


Point::Point(double x_ /*= 0*/, double y_ /*= 0*/)
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


Basis::Basis(size_t length_ /*= 0*/)
{
    length = length_;
    polynomials = new Polynomial[length];
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
