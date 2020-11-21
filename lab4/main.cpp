#include "Header.hpp"

int main()
{
    size_t n = 4;
    int size = n + 1;
    double leftBorder = -1;
    double rightBorder = 1;

    Grid uniformGrid = Grid(size);
    Grid ChebishevGrid = Grid(size);
    
    Polynomial p1 = Polynomial(size);
    Polynomial pBuffer = Polynomial(size);
    
    Basis basis = Basis(n);


    MakeMesh(leftBorder, rightBorder, uniformGrid, 1);
    MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2);

    ReadInit();
    AllocateMemory();
    ReadData();

    LagrangeInterpolate(uniformGrid, basis, p1, pBuffer);
    SplineInterpolate();

    WriteData();
    FreeMemory();

    std::cout << "Hello world xdDDDddD\n";

    return 0;
}
