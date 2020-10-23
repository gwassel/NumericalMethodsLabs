#include <iostream>
#include "FuncsDefs.h"
#include "header.hpp"
#define my_type double
//#define my_type float	
int main()
{
    std::cout << "------------------------------------------------------Molochkov------------------------------------------------------" << std::endl;
    size_t n = 4;

    //input data
    double** matrixA = nullptr;
    double* vectorB = nullptr;

    double thau = 0;
    double epsilon = 0;

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

    std::cout << "--------------------------------------------------------Emelin-------------------------------------------------------" << std::endl;
    const int Esize = 204;

    my_type** EmatrixA;
    my_type** EcopyA;
    my_type** EDplusL;
    my_type** EDplusLI;
    my_type** EBI;
    my_type** EmBuffer;
    my_type** ET;
    my_type** EQ;
    my_type** ER;

    my_type* EmatrixC_L;//нижняя
    my_type* EmatrixC_D;//диагональ
    my_type* EmatrixC_U;//верхняя

    my_type* EvectorB;
    my_type* EvectorBCopy;
    my_type* EBF;
    my_type* EvectorD;
    my_type* EvectorX1;
    my_type* EvectorX2;
    my_type* EvectorXZ;
    my_type* EvectorXR;
    my_type* EcopyB;
    my_type* Ebuffy;

    const std::string EpathData = "data/";

    std::string EfileNameMatrix;
    std::string EfileNameVector;
    std::string EfolderName;

    my_type EnormC_L;
    my_type EnormC_D;
    my_type EnormC_U;

    EfolderName = "results/res" + std::to_string(2);
    EfileNameMatrix = EpathData + "Ematrix3.txt";
    EfileNameVector = EpathData + "Evector3.txt";

    const std::string EfileMatrixAName = "/matrixA";
    const std::string EfileVectorXZName = "/vectorXZ";
    const std::string EfileVectorXRName = "/vectorXR";
    const std::string EfileVectorBName = "/vectorB";
    EWriteMatrixAB(Esize);

    EAllocateMemory(EmatrixA, EmatrixC_L, EmatrixC_D, EmatrixC_U, EvectorB, EvectorBCopy, EvectorX1, EvectorX2, EvectorXZ,
                    EvectorXR,EvectorD, EBF, EcopyA, EDplusL, EDplusLI,EBI,EmBuffer,ET,EQ,ER, EcopyB, Ebuffy, Esize);

    EReadData(EfileNameMatrix, EfileNameVector,EmatrixA,EvectorB, EvectorBCopy, EcopyA, EDplusL, EcopyB, Esize);

    ECalculations(EmatrixA, EvectorB, EvectorBCopy, EvectorX1, EvectorX2, EvectorXZ, EvectorXR, EmatrixC_U, EmatrixC_D,
                  EmatrixC_L, EvectorD, EBF, EnormC_L, EnormC_D, EnormC_U, EcopyA, EDplusL, EDplusLI, EBI, EmBuffer, ET, EQ, ER, EcopyB, Ebuffy, Esize);

    EWriteData(EfolderName + EfileMatrixAName, EfolderName + EfileVectorBName,
              EfolderName + EfileVectorXZName, EfolderName + EfileVectorXRName, EmatrixA, EvectorXZ, EvectorXR, EvectorB,Esize);
    EFreeMemory(EmatrixA, EmatrixC_L, EmatrixC_D, EmatrixC_U, EvectorB, EvectorBCopy, EvectorX1, EvectorX2, EvectorXZ,
        EvectorXR, EvectorD, EBF, EcopyA, EDplusL, EDplusLI, EBI, EmBuffer, ET, EQ, ER, EcopyB, Ebuffy, Esize);







}
