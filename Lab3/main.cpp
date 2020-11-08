#include "Header.hpp"
#include "../nlohmann/json.hpp"

//SEV - searching eigen values
int main()
{
    const int n = 4;
    const int Esize = 4;
    const double accuracy = 1e-5;

    double** matrixA;
//E vars
    double** EmatrixT;
    double** EmatrixQ;
    double** EmatrixR;
    double** EmatrixC;
    double** EmatrixBuffer1;
    double** EmatrixBuffer2;
    double** EmatrixEigenVectors;

    double* Elambda;
    double* EvectorX1;
    double* EvectorX2;
    double* EvectorBstar;
    double* EvectorBuffer;
    double* EvectorLambda;
//M vars
    double** MmatrixH;
    double** MmatrixAk;
    double** MmatrixQ;
    double** MmatrixR;
    double** MmatrixBuffer;
    double* MvectorLambdaOld;
    double* MvectorLambdaNew;

    nlohmann::json Ej;
    std::ifstream Efile("menu/Econfig.json");
    const int ESYSTEM = 0;
    Efile >> Ej;

    const std::string EfileNameMatrix = Ej[ESYSTEM]["package_init"].get<std::string>() + Ej[ESYSTEM]["matrix_init_name"].get<std::string>();
    const std::string EfileNameEigenValsInit = Ej[ESYSTEM]["package_init"].get<std::string>() + Ej[ESYSTEM]["eigen_values_file_init_name"].get<std::string>();

    const std::string EfileMatrixEVecName = Ej[ESYSTEM]["package_res"].get<std::string>() + Ej[ESYSTEM]["eigen_vectors_file_res_name"].get<std::string>();
    const std::string EfileMatrixEValName = Ej[ESYSTEM]["package_res"].get<std::string>() + Ej[ESYSTEM]["eigen_values_file_res_name"].get<std::string>();


    EAllocateMemory(matrixA,EvectorX1,EvectorX2,EvectorLambda,EvectorBstar,EvectorBuffer,EmatrixBuffer1,
                    EmatrixBuffer2,EmatrixT,EmatrixQ,EmatrixR,EmatrixEigenVectors,EmatrixC,Esize);
    MAllocateMemory(MmatrixH, MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, n);


    ReadData(EfileNameMatrix, EfileNameEigenValsInit, matrixA, EvectorLambda, Esize);

    HessenbergForm(matrixA, MmatrixH, n);
    
    PrintMatrix(matrixA, n, "matrix A");
    PrintMatrix(MmatrixH, n, "matrixH");

    MatrixCopy(MmatrixAk, matrixA, n);
    int k1 = SimpleQRIterations(MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, accuracy, n);
    std::cout << "iter for A: " << k1 << "\n";
    PrintVector(MvectorLambdaNew, n, "vectorLambda for matrixA");
    int k2 = SimpleQRIterations(MmatrixH, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, accuracy, n);
    std::cout << "iter for H: " << k2 << "\n";
    PrintVector(MvectorLambdaNew, n, "vectorLambda for matrixH");



    std::cout << "--------------------------------------------------------Emelin-------------------------------------------------------" << std::endl;

    ECalculations(matrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);

    EWriteData(EfileMatrixEVecName, EfileMatrixEValName, EmatrixEigenVectors, EvectorLambda, Esize);

    EFreeMemory(matrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);
    MFreeMemory();

    system("pause");
}
