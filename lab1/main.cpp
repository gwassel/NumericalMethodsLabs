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
    my_type** matrixAInverted = nullptr;

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
    std::string fileNameAInv = "";
    std::string fileNameA_AInv = "";
    
    size_t n = 0;

    ReadInit(pathConfig, fileNameA, fileNameB, fileNameQ, fileNameR, fileNameX, fileNameAInv, fileNameA_AInv, n);
    
    AllocateMemory(matrixA, matrixT, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1,
            matrixBuffer2, matrixAInverted, vectorBuffer, n);

    ReadData(fileNameA, fileNameB, matrixA, vectorB, n);   
   
    QRCalculations(matrixA, matrixT, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1, matrixBuffer2, vectorBuffer, n); 

    MatrixInverseTR(matrixT, matrixR, matrixAInverted, matrixBuffer1, n);

    MatrixMult(matrixA, matrixAInverted, matrixBuffer1, n);
    MatrixMult(matrixA, vectorX, vectorBuffer, n);// vectorBuffer=B*

    for(int i = 0; i < n; ++i)
    {
        vectorBuffer[i] = vectorBuffer[i] - vectorB[i];
    }
    std::cout << "octaedral vector norm: " << OctahedralVectorNorm(vectorBuffer, n) << std::endl;
    std::cout << "cubic vector norm" << CubicVectorNorm(vectorBuffer, n) << std::endl;

    std::cout << "cond octaedral " << OctahedralMatrixNorm(matrixA, n) * OctahedralMatrixNorm(matrixAInverted, n) << std::endl;
    std::cout << "cond cubic " << CubicMatrixNorm(matrixA, n) * CubicMatrixNorm(matrixAInverted, n) << std::endl;

    WriteData(fileNameQ, fileNameR, fileNameX, fileNameA_AInv, fileNameAInv, matrixQ, matrixR, vectorX, matrixBuffer1, matrixAInverted, n);
    WriteVector("data/vectorBStarred", "BStarred", vectorBuffer, n);
    FreeMemory(matrixA, matrixT, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1,
            matrixBuffer2, matrixAInverted, vectorBuffer, n);

    return 0;
}
