#include "funcs.hpp"

//is correct funcs
//is Q*transp(Q) ==1
//is QR=A
//is R a upper triangle

//init

int main()
{
    //QR variables
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
    
    std::string fileNameA = "";
    std::string fileNameB = "";
    
    std::string fileNameQ = "";
    std::string fileNameR = "";
    std::string fileNameX = "";
    std::string fileNameAInv = "";
    std::string fileNameA_AInv = "";
    
    size_t n = 0;

    //Gauss variables
    const int Esize = 4;
    my_type** EA;
    my_type** EIA;
    my_type* EB;
    my_type* EBF;
    my_type* EdB;
    my_type* EX1;
    my_type* EX2;
    my_type Econd;
    my_type Edx;
    my_type Edb;
    bool EnoProblem;
    std::tuple<my_type, my_type, my_type> Ed11;
    std::tuple<my_type, my_type> Ec11;

    const std::string EpathConfig = "configs/config.json";
    const std::string EpathData = "data/";
    std::string EfileNameMatrix;
    std::string EfileNameVector;
    const std::string EfileMatrixAName = "/matrixA";
    const std::string EfileMatrixAIName = "/matrixAI";
    const std::string EfileVectorX1Name = "/vectorX1";
    const std::string EfileVectorX2Name = "/vectorX2";
    const std::string EfileVectorBName = "/vectorB";
    const std::string EfileParamsName = "/Params";
    std::string EfolderName;

    ReadInit(pathConfig, fileNameA, fileNameB, fileNameQ, fileNameR, fileNameX, fileNameAInv, fileNameA_AInv, n);

    MMain(matrixA, vectorB, vectorX, matrixQ, matrixR, matrixT, matrixAInverted, matrixBuffer1, matrixBuffer2, vectorBuffer,
            fileNameA, fileNameB, fileNameQ, fileNameR, fileNameX, fileNameAInv, fileNameA_AInv, n);
    
    EMain(Esize, EA, EIA, EB, EBF, EdB, EX1, EX2, Econd, Edx, Edb, EnoProblem, Ed11, Ec11, EpathConfig, EpathData,
            EfileNameMatrix, EfileNameVector, EfileMatrixAName, EfileMatrixAIName, EfileVectorX1Name, EfileVectorX2Name,
            EfileVectorBName, EfileParamsName, EfolderName);

    return 0;
}
