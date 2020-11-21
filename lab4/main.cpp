#include "Header.hpp"

int main()
{
    size_t n = 3;

    Grid uniformGrid = Grid(n);
    Grid ChebishevGrid = Grid(n);
    
    Polynomial p1 = Polynomial(n);
    Polynomial pBuffer = Polynomial(n);
    
    Basis basis = Basis(n);


    ReadInit();
    AllocateMemory();
    ReadData();

    MakeUniformMesh();
    MakeChebishevMesh();//name may be incorrect

    LagrangeInterpolate(uniformGrid, basis, p1, pBuffer);
    SplineInterpolate();

    WriteData();
    FreeMemory();

    std::cout << "Hello world xdDDDddD\n";

    return 0;
}
