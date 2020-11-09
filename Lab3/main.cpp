#include "Header.hpp"
#include "../nlohmann/json.hpp"

//SEV - searching eigen values
int main()
{
    const int n = 4;
    const int Esize = n;
    const double accuracy = 1e-6;

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


//readinit
    nlohmann::json Ej;
    std::ifstream Efile("menu/Econfig.json");
    const int ESYSTEM = 0;
    Efile >> Ej;

    const std::string EfileNameMatrix = Ej[ESYSTEM]["package_init"].get<std::string>() + Ej[ESYSTEM]["matrix_init_name"].get<std::string>();
    const std::string EfileNameEigenValsInit = Ej[ESYSTEM]["package_init"].get<std::string>() + Ej[ESYSTEM]["eigen_values_file_init_name"].get<std::string>();

    const std::string EfileMatrixEVecName = Ej[ESYSTEM]["package_res"].get<std::string>() + Ej[ESYSTEM]["eigen_vectors_file_res_name"].get<std::string>();
    const std::string EfileMatrixEValName = Ej[ESYSTEM]["package_res"].get<std::string>() + Ej[ESYSTEM]["eigen_values_file_res_name"].get<std::string>();
//readinit


    EAllocateMemory(matrixA,EvectorX1,EvectorX2,EvectorLambda,EvectorBstar,EvectorBuffer,EmatrixBuffer1,
                    EmatrixBuffer2,EmatrixT,EmatrixQ,EmatrixR,EmatrixEigenVectors,EmatrixC,Esize);
    MAllocateMemory(MmatrixH, MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, n);


    ReadData(EfileNameMatrix, EfileNameEigenValsInit, matrixA, EvectorLambda, Esize);


//Mcalc
    HessenbergForm(matrixA, MmatrixH, n);
    
    PrintMatrix(matrixA, n, "matrix A");
    PrintMatrix(MmatrixH, n, "matrixH");

    MatrixCopy(MmatrixAk, matrixA, n);
    int k1 = SimpleQRIterations(MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, accuracy, n);
    std::cout << "iter for A: " << k1 << "\n";
    PrintVector(MvectorLambdaNew, n, "vectorLambdaNew for matrixA");

    MatrixCopy(MmatrixAk, MmatrixH, n);
    int k2 = SimpleQRIterations(MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, accuracy, n);
    std::cout << "iter for H: " << k2 << "\n";
    PrintVector(MvectorLambdaNew, n, "vectorLambda for matrixH");

    MatrixCopy(MmatrixAk, matrixA, n);
    int k3 = ShiftQRIterations(MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, accuracy, n);
    PrintVector(MvectorLambdaNew, n, "vectorLambdaNew for matrixA with shafts");
    PrintVector(MvectorLambdaOld, n-1, "number of iterations by steps");
    
    MatrixCopy(MmatrixAk, MmatrixH, n);
    int k4 = ShiftQRIterations(MmatrixAk, MmatrixQ, MmatrixR, MmatrixBuffer, MvectorLambdaOld, MvectorLambdaNew, accuracy, n);
    PrintVector(MvectorLambdaNew, n, "vectorLambda for matrixH with shafts");
    PrintVector(MvectorLambdaOld, n-1, "number of iterations by steps");
//Mcalc
    
    // EvectorLambda[0] = MvectorLambdaNew[1];
    // EvectorLambda[1] = MvectorLambdaNew[2];
    // EvectorLambda[2] = MvectorLambdaNew[0];
    // EvectorLambda[3] = MvectorLambdaNew[3];

    for(int i = 0; i < n; ++i)
        EvectorLambda[i]=MvectorLambdaNew[i];

    std::cout << "--------------------------------------------------------Emelin-------------------------------------------------------" << std::endl;

    ECalculations(matrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);

    EWriteData(EfileMatrixEVecName, EfileMatrixEValName, EmatrixEigenVectors, EvectorLambda, Esize);

    MatrixMultVector(matrixA, EmatrixEigenVectors[3], EvectorBuffer, n);
    for(int i = 0; i < n; ++i)
        EmatrixEigenVectors[3][i] *= EvectorLambda[3];

    PrintVector(EvectorBuffer, n, "this need to be eq to that");
    PrintVector(EmatrixEigenVectors[0], n, "this need to be eq to that");
    

    EFreeMemory(matrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);
    MFreeMemory();

    system("pause");
}
