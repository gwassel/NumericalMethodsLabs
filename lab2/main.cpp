#include "header.hpp"

int main()
{
    int n = 4;

    //input data
    double** matrixA = nullptr;
    double* vectorB = nullptr;

    //input files names
    const std::string FILE_NAME_CFG = "cfg/config.json";
    std::string fileNameMatrixA = "data2/matrix.dat";
    std::string fileNameVectorB = "data2/vector.dat";
    std::string fileNameVectorX0 = "data2/init_approx";

    std::string fileNameMatrixC = "results/C";
    std::string fileNameVectorX = "results/X";
    std::string fileNameVectorY = "results/Y";

    //MData
    double** matrixC = nullptr;
    double* vectorXCurrent = nullptr;
    double* vectorXFollow = nullptr;
    double* vectorY = nullptr;

    //MBuffer
    double* vectorBuffer = nullptr;


//    ReadInit();

    MAllocateMemory(matrixA, vectorB, matrixC, vectorXCurrent, vectorXFollow, vectorY, n);
    MAllocateBuffer(vectorBuffer, n);

    ReadData(fileNameMatrixA, fileNameVectorB, fileNameVectorX0, matrixA, vectorB, vectorXCurrent, n);
    
    MCalculations(matrixA, vectorB, matrixC, vectorXCurrent, vectorXFollow, vectorY, vectorBuffer, n);

    MWriteData(fileNameMatrixC, fileNameVectorX, fileNameVectorY, matrixC, vectorXCurrent, vectorY, n);

    MFreeMemory(matrixA, vectorB, matrixC, vectorXCurrent, vectorXFollow, vectorY, n);
    MFreeBuffer(vectorBuffer, n);

    return 0;
}
