#include "header.hpp"

bool IfStop(const double* const* matrixC, const double* vectorXCurrent, const double* vectorXFollow, double* &vectorBuffer, const double epsilon, const size_t n)
{
    // std::cout << "IfStop\n";
    VectorDiff(vectorXCurrent, vectorXFollow, vectorBuffer, n);
    double normDeltaX = CubicVectorNorm(vectorBuffer, n);
    double normC = CubicMatrixNorm(matrixC, n);

    // std::cout << "normdX = " << normDeltaX << " normC = " << normC << "\n";

    bool ifStop = (normDeltaX <= (1 - normC) / normC * epsilon);
    // std::cout << ifStop << "\n";
    return ifStop;
}

void GetCY_SimpleIter(const double* const* matrixA, const double* vectorB, const double thau, double** &matrixC, double* &vectorY, const size_t n)
{
    std::cout << "GetCY simple iter\n";
    MatrixCopy(matrixC, matrixA, n);
    MatrixDot(matrixC, thau, n);
    MatrixEDiff(matrixC, n);
    VectorCopy(vectorY, vectorB, n);
    VectorDot(vectorY, thau, n);
    std::cout << "norm C: " << CubicMatrixNorm(matrixC, n) << "\n";
}

double ResidualNorm(const double* const* matrixA, const double* vectorB, const double* vectorX, double* &vectorBuffer1, double* &vectorBuffer2, const size_t n)
{
    MatrixMult(matrixA, vectorX, vectorBuffer1, n);
    VectorDiff(vectorBuffer1, vectorB, vectorBuffer2, n);
    
    return CubicVectorNorm(vectorBuffer2, n);
}

void GetCY_Jacobi(const double* const* matrixA, const double* vectorB, double** &matrixC, double* &vectorY, const size_t n)
{
    std::cout << "GetCY Jacobi\n";
    for(int i = 0; i < n; ++i)
    {
        double Aii = matrixA[i][i];
        for(int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                matrixC[i][j] = 0;
            }
            else
            {
                matrixC[i][j] = -matrixA[i][j] / Aii;
            }
        }
        vectorY[i] = vectorB[i] / Aii;
    }
    std::cout << "norm C: " << CubicMatrixNorm(matrixC, n) << "\n";
}

void Iterations(const double* const* matrixC, double* &vectorXCurrent, double* &vectorXFollow, const double* vectorY, double* &vectorBuffer, const double epsilon, const size_t n)
{
    VectorCopy(vectorXFollow, vectorXCurrent, n);
    std::cout << "Iterations\n";
    
    int iterationNumber = 1;

    do
    {
//        std::cout << "\niteration number = " << iterationNumber << "\n\n";
        iterationNumber++;

        VectorCopy(vectorXCurrent, vectorXFollow, n); // Xk=Xk+1
        MatrixMult(matrixC, vectorXCurrent, vectorBuffer, n); // Cx
//        LogVector("Cx:", vectorBuffer, n); 


        VectorAdd(vectorBuffer, vectorY, vectorXFollow, n); // Xk+1 = Cxk + y

//        LogVector("X:", vectorXCurrent, n);
//        LogVector("Xnext:", vectorXFollow, n);
    }
    while(!IfStop(matrixC, vectorXCurrent, vectorXFollow, vectorBuffer, epsilon, n) && iterationNumber < 100000);
    LogVector("x:", vectorXCurrent, n);
    std::cout << "iterations: " << iterationNumber << "\n";
}

void MCalculations(const double* const* matrixA, const double* vectorB, double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorY, double* &vectorBuffer, const double thau, const double epsilon, const size_t n)
{
    std::cout << "MCalculations\n";
    //GetCY_SimpleIter(matrixA, vectorB, thau, matrixC, vectorY, n);
    GetCY_Jacobi(matrixA, vectorB, matrixC, vectorY, n);
    Iterations(matrixC, vectorXCurrent, vectorXFollow, vectorY, vectorBuffer, epsilon, n);
    std::cout << "res norm: " << ResidualNorm(matrixA, vectorB, vectorXCurrent, vectorBuffer, vectorXFollow, n) << "\n";
}

