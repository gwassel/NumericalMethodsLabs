#include "Header.hpp"
#include "../nlohmann/json.hpp"

//SEV - searching eigen values
int main()
{
    const int n = 4;
    const int Esize = 4;

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
    double** matrixH;


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

    //MAllocateMemory(matrixH, n);
    AllocateMatrix(matrixH, n);
    ReadData(EfileNameMatrix, EfileNameEigenValsInit, matrixA, EvectorLambda, Esize);

    HessenbergForm(matrixA, matrixH, n);
    
    PrintMatrix(A, n, "matrix A");
    PrintMatrix(H, n, "matrixH");

    


    std::cout << "--------------------------------------------------------Emelin-------------------------------------------------------" << std::endl;

    ECalculations(matrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);

    EWriteData(EfileMatrixEVecName, EfileMatrixEValName, EmatrixEigenVectors, EvectorLambda, Esize);

    EFreeMemory(matrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);

    system("pause");
}
