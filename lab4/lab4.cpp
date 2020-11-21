#include <iostream>
#include "Header.h"


int main()
{
    int n = 4;
    int size = n + 1;
    double* PointVector = new double[size];
    double A = -1;
    double B = 1;

    MakeMesh(A, B, PointVector, size, 2);

    for (int i = 0; i < size; i++) {
        std::cout << PointVector[i] << std::endl;
    }

    //declare variables
/*
    ReadInit();
    AllocateMemory();
    ReadData();

    MakeUniformMesh();
    MakeChebishevMesh();//name may be incorrect

    LagrangeInterpolate();
    SplineInterpolate();

    WriteData();
    FreeMemory();*/
    return 0;
}
