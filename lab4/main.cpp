#include <iostream>

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
    return 0;
}
