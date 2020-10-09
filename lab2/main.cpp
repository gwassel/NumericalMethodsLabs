#include "header.hpp"

int main()
{
    size_t n = 4;

    //input data
    double** matrixA = nullptr;
    double* vectorB = nullptr;
    
    double thau = 0.1;
    double epsilon = 0.0001;
    
    //input files names
    const std::string FILE_NAME_CFG = "cfg/config.json";

    std::string folderInput = "data/";
    std::string fileNameMatrixA = "matrix.dat";
    std::string fileNameVectorB = "vector.dat";
    std::string fileNameVectorX0 = "init_approx";

    std::string folderOutput = "results/";
    
    //test
    std::string folderTest = "test/";
    std::string testA = "";
    std::string testB = "";
    std::string testSolve = "";
    std::string testInitApprox = "";

    double** matrixAT = nullptr;
    double* vectorBT = nullptr;
    double* vectorSolveT = nullptr;
    double* vectorSolveApprox = nullptr;
    

    //MData
    double** matrixC = nullptr;
    double* vectorXCurrent = nullptr;
    double* vectorXFollow = nullptr;
    double* vectorY = nullptr;

    //MBuffer
    double* vectorBuffer = nullptr;


    ReadInit(FILE_NAME_CFG, folderInput, fileNameMatrixA, fileNameVectorB, fileNameVectorX0, folderOutput, folderTest, testA, testB, testSolve, testInitApprox, thau, epsilon, n);

    MAllocateMemory(matrixA, vectorB, matrixC, vectorXCurrent, vectorXFollow, vectorY, n);
    MAllocateBuffer(vectorBuffer, n);

    ReadData(folderInput, fileNameMatrixA, fileNameVectorB, fileNameVectorX0, matrixA, vectorB, vectorXCurrent, n);
    
    MCalculations(matrixA, vectorB, matrixC, vectorXCurrent, vectorXFollow, vectorY, vectorBuffer, thau, epsilon, n);

    MWriteData(folderOutput, matrixC, vectorXCurrent, vectorY, n);

    MFreeMemory(matrixA, vectorB, matrixC, vectorXCurrent, vectorXFollow, vectorY, n);
    MFreeBuffer(vectorBuffer, n);

    return 0;
}
