#include "funcs.hpp"

#define my_type double


//is correct funcs
//is Q*transp(Q) ==1
//is QR=A
//is R a upper triangle

//init

int main()
{
    my_type** matrixA = nullptr;
    my_type* vectorB = nullptr;
    my_type* vectorX = nullptr;
    my_type** matrixQ = nullptr;
    my_type** matrixR = nullptr;
    my_type** matrixT = nullptr;
    
    my_type** matrixBuffer1 = nullptr;
    my_type** matrixBuffer2 = nullptr;
    my_type* vectorBuffer = nullptr;

    const std::string pathConfig = "configs/config.json";
    const std::string pathData = "data/";
    
    //inputs
    std::string fileNameA = "";
    std::string fileNameB = "";
    
    //outputs
    std::string fileNameQ = "";
    std::string fileNameR = "";
    std::string fileNameX = "";
    
    size_t n = 0;

    ReadInit(pathConfig, fileNameA, fileNameB, fileNameQ, fileNameR, fileNameX, n);
    
    AllocateMemory(matrixA, matrixT, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1,
            matrixBuffer2, vectorBuffer, n);

    ReadData(fileNameA, fileNameB, matrixA, vectorB, n);   
   
    QRCalculations(matrixA, matrixT, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1, matrixBuffer2, vectorBuffer, n); 

    MatrixInverseTR(matrixT, matrixR, matrixBuffer1, matrixBuffer2, n);

    WriteData(fileNameQ, fileNameR, fileNameX, matrixQ, matrixR, vectorX, n);

    FreeMemory(matrixA, matrixT, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1,
            matrixBuffer2, vectorBuffer, n);

    return 0;
}
