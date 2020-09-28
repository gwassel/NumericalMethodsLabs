#include "header.hpp"

bool IfStop(double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorBuffer, const double epsilon, const size_t n)
{
    std::cout << "IfStop\n";
    VectorDiff(vectorXCurrent, vectorXFollow, vectorBuffer, n);
    double normDeltaX = CubicVectorNorm(vectorBuffer, n);
    double normC = CubicMatrixNorm(matrixC, n);

    std::cout << "normdX = " << normDeltaX << " normC = " << normC << "\n";

    bool ifStop = (normDeltaX <= (1 - normC) / normC * epsilon);
    std::cout << ifStop << "\n";
    return ifStop;
}

void GetCY(double** &matrixA, double* &vectorB, double thau, double** &matrixC, double* &vectorY, const size_t n)
{
    std::cout << "GetCY\n";
    MatrixCopy(matrixC, matrixA, n);
    MatrixDot(matrixC, thau, n);
    MatrixEDiff(matrixC, n);
    VectorCopy(vectorY, vectorB, n);
    VectorDot(vectorY, thau, n);
}

void Iterations(double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorY, double* &vectorBuffer, const double epsilon, const size_t n)
{
    VectorCopy(vectorXFollow, vectorXCurrent, n);
    std::cout << "Iterations\n";
    
    int iterationNumber = 1;

    do
    {
        std::cout << "\niteration number = " << iterationNumber++ << "\n\n";

        VectorCopy(vectorXCurrent, vectorXFollow, n); // Xk=Xk+1
        MatrixMult(matrixC, vectorXCurrent, vectorBuffer, n); // Cx
        LogVector("Cx:", vectorBuffer, n); 


        VectorAdd(vectorBuffer, vectorY, vectorXFollow, n); // Xk+1 = Cxk + y

        LogVector("X:", vectorXCurrent, n);
        LogVector("Xnext:", vectorXFollow, n);
    }
    while(!IfStop(matrixC, vectorXCurrent, vectorXFollow, vectorBuffer, epsilon, n) && iterationNumber < 1000);
}

void MCalculations(double** &matrixA, double* &vectorB, double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorY, double* &vectorBuffer, const size_t n)
{
    std::cout << "MCalculations\n";
    GetCY(matrixA, vectorB, 0.00000001, matrixC, vectorY, n);
    Iterations(matrixC, vectorXCurrent, vectorXFollow, vectorY, vectorBuffer, 0.1, n);
}



