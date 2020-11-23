#include "Header.hpp"

int main()
{
    //init variables
    size_t n = 4;
    size_t size = n + 1;
    size_t testSize = 1024 + 1;

    double leftBorder = -1;
    double rightBorder = 1;

    //non init variables
    Grid uniformGrid = Grid(size);
    Grid ChebishevGrid = Grid(size);
    Grid testGrid = Grid(testSize);

    Polynomial p1 = Polynomial(size);
    Polynomial p2 = Polynomial(size);
    Polynomial pBuffer = Polynomial(size);
    Polynomial monomial = Polynomial(2);
    
    Basis basis1 = Basis(size);
    Basis basis2 = Basis(size);

    TridiagonalMatrix splineMatrix;
    double (*f)(double) = f2; //init function to interpolate

    MakeMesh(leftBorder, rightBorder, uniformGrid, 1, f);
    MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2, f);

    ReadInit();
    ReadData();

    LagrangeInterpolate(uniformGrid, basis1, p1, pBuffer, monomial, "Lagrange Polynomial on uniform grid");
    LagrangeInterpolate(ChebishevGrid, basis2, p2, pBuffer, monomial, "Lagrange Polynomial on Chebishev grid");

    SplineInterpolate();

    WriteData();

    test(p1, uniformGrid, "Lagrange on uniform grid");
    test(p2, ChebishevGrid, "Lagrange on Chebishev grid");
    double err = CountError(p1, testGrid, f, leftBorder, rightBorder);
    std::cout << "err on " << testSize - 1 << " grid = " << err << "\n";
    return 0;
}
