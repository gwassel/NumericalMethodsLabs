#include <cmath>
#include <utility>
#include <iostream>
#include "Header.hpp"

void MAllocateMemory(double** &matrixH, double** &matrixAk, double** &matrixQ, double** &matrixR, double** &matrixBuffer, double* &vectorLambdaOld, 
        double* &vectorLambdaNew, const size_t n)
{
    AllocateMatrix(matrixH, n);
    AllocateMatrix(matrixAk, n);
    AllocateMatrix(matrixQ, n);
    AllocateMatrix(matrixR, n);
    AllocateMatrix(matrixBuffer, n);

    AllocateVector(vectorLambdaOld, n);
    AllocateVector(vectorLambdaNew, n);
}

void MFreeMemory()
{
    std::cout << "Need to add MFreeMemory!!\n";
}

void HessenbergForm(double** &matrixA, double** &matrixH, const size_t n)
{
    MatrixCopy(matrixH, matrixA, n);
    double c, s, radical, c1, c2;

    for(int j = 0; j < n - 2; ++j)
    {
        for(int i = j + 2; i < n; ++i)
        {
            c = matrixH[j+1][j];
            s = matrixH[i][j];
            
            radical = 1 / sqrt(c*c + s*s);

            c *= radical;
            s *= radical;

            for(int k = j; k < n; ++k)
            {
                c1 = c * matrixH[j+1][k] + s * matrixH[i][k];
                c2 = (-s) * matrixH[j+1][k] + c * matrixH[i][k];

                matrixH[j+1][k] = c1;
                matrixH[i][k] = c2;
            }

            for(int k = 0; k < n; ++k)
            {
                c1 = c * matrixH[k][j+1] + s * matrixH[k][i];
                c2 = (-s) * matrixH[k][j+1] + c * matrixH[k][i];

                matrixH[k][j+1] = c1;
                matrixH[k][i] = c2;
            }
        }
    }
}

int SimpleQRIterations(double** &matrixAk, double** &matrixQ, double** &matrixR, double** &matrixBuffer, double* &vectorLambdaOld, 
        double* &vectorLambdaNew, double accuracy, const size_t n)//Ak=A0
{
    int k = 0;
    double delta = 100;
    PrintMatrix(matrixAk, n, "matrixAk in SimpleQRI");

    while(k < 100 && delta > accuracy)
    {
        QRDecomposerLite(matrixAk, matrixBuffer, matrixQ, matrixR, n);
        MatrixMultMatrix(matrixR, matrixQ, matrixAk, n);

        for(int i = 0; i < n; ++i){
            vectorLambdaOld[i] = vectorLambdaNew[i];
            vectorLambdaNew[i] = matrixAk[i][i];
            matrixBuffer[0][i] = vectorLambdaNew[i] - vectorLambdaOld[i];
        }

        delta = CubicVectorNorm(matrixBuffer[0], n);
        ++k;
    }

    return k;
}
