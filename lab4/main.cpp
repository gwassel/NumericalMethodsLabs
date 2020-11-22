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
    Polynomial p2 = Polynomial(size);
    Polynomial pBuffer = Polynomial(size);
    Polynomial monomial = Polynomial(2);
    
    Basis basis1 = Basis(size);
    Basis basis2 = Basis(size);

    TridiagonalMatrix splineMatrix;


    MakeMesh(leftBorder, rightBorder, uniformGrid, 1, f2);
    MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2, f2);

    ReadInit();
    AllocateMemory();
    ReadData();

    LagrangeInterpolate(uniformGrid, basis1, p1, pBuffer, monomial, "Lagrange Polynomial on uniform grid");
    LagrangeInterpolate(ChebishevGrid, basis2, p2, pBuffer, monomial, "Lagrange Polynomial on Chebishev grid");

    SplineInterpolate();

    WriteData();
    FreeMemory();

    test(p1, uniformGrid, "Lagrange on uniform grid");
    test(p2, ChebishevGrid, "Lagrange on Chebishev grid");

    return 0;
}
