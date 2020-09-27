#include "header.hpp"

int main()
{
    //input data
    double** matrixA = nullptr;
    double* vectorB = nullptr;

    //input files names
    std::string matrixFile = "";
    std::string vectorFile = "";


    ReadInit();

    MAllocateMemory();
    MAllocateBuffer();

    MCalculations();

    MWriteData();

    MFreeMemory();
    MFreeBuffer();
    return 0;
}
