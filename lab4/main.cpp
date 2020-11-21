#include "Header.hpp"

int main()
{
    //declare variables

    ReadInit();
    AllocateMemory();
    ReadData();

    MakeUniformMesh();
    MakeChebishevMesh();//name may be incorrect

    LagrangeInterpolate();
    SplineInterpolate();

    WriteData();
    FreeMemory();

    std::cout << "Hello world xdDDDddD\n";

    return 0;
}
