#include "Header.hpp"

int main()
{
    size_t n = 4;
    int size = n + 1;
    double* PointVector = new double[size];
    double leftBorder = -1;
    double rightBorder = 1;

    Grid uniformGrid = Grid(size);
    Grid ChebishevGrid = Grid(size);
    
    Polynomial p1 = Polynomial(n);
    Polynomial pBuffer = Polynomial(n);
    
    Basis basis = Basis(n);


    MakeMesh(leftBorder, rightBorder, uniformGrid, 1);
    MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2);

    for (int i = 0; i < size; i++) {
        std::cout << ChebishevGrid.points[i].x << std::endl;
    }


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
