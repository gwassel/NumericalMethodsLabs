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

    Polynomial monomial = Polynomial(2);
    
    Basis basis = Basis(size);


    MakeMesh(leftBorder, rightBorder, uniformGrid, 1);
    //MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2);

    ReadInit();
    AllocateMemory();
    ReadData();

    LagrangeInterpolate(uniformGrid, basis, p1, pBuffer, monomial);

    double err = test(p1, uniformGrid);
    SplineInterpolate();

    for(int i = 0; i < p1.length; ++i)
        std::cout << p1.coefficents[i] << "x^" << i << " ";
    std::cout << "\n";

    WriteData();
    FreeMemory();

    std::cout << "err:" << err << "\n";

    return 0;
}
