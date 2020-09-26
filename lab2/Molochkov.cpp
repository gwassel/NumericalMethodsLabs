#include "header.hpp"

bool IfStop(double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorBuffer, const double epsilon, const size_t n)
{
    VectorDiff(vectorXCurrent, vectorXFollow, vectorBuffer, n);
    double normDeltaX = CubicVectorNorm(vectorBuffer, n);
    double normC = CubicMatrixNorm(matrixC, n);
    return (normDeltaX <= (1 - normC) / normC * epsilon);
}

void Iterations(double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorY, double* &vectorBuffer, const double epsilon, const size_t n)
{
    do
    {
        MatrixMult(matrixC, vectorXCurrent, vectorBuffer, n); // Cx
        VectorCopy(vectorXCurrent, vectorXFollow, n); // Xk=Xk+1
        VectorAdd(vectorBuffer, vectorY, vectorXFollow, n); // Xk+1 = Cxk + y
    }
    while(ifStop(matrixC, vectorXCurrent, vectorXFollow, vectorBuffer, epsilon, n))
}

void GetCY(double** &matrixA, double* &vectorB, double thau, double** &matrixC, double* &vectorY, const size_t n)
{
    MatrixEDiff(matrixA, matrixC, n);
    MatrixCopy(matrixY, matrixB, n);
}
