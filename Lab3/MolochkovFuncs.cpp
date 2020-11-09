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

    while(k < 2000 && delta > accuracy)
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

int ShiftQRIterations(double** &matrixAk, double** &matrixQ, double** &matrixR, double** &matrixBuffer, double* &vectorOfK, 
        double* &vectorLambdaNew, double accuracy, const size_t n)
{
    int k = 0;
    double sigma = 0, delta = 100;
    for(int i = n; i > 1; --i){
        k = 0;
        delta = 500;
        while(k < 5000 && delta > accuracy)
        {
            sigma = matrixAk[i-1][i-1];
            MatrixResE(matrixAk, i, sigma);
            QRDecomposerLite(matrixAk, matrixBuffer, matrixQ, matrixR, i);
            MatrixMultMatrix(matrixR, matrixQ, matrixAk, i);
            MatrixResE(matrixAk, i, -sigma);
            
            for(int j = 0; j < i - 1; ++j)
                matrixBuffer[0][j] = matrixAk[i-1][j];

            delta = CubicVectorNorm(matrixBuffer[0], i);
            ++k;
        }
        vectorOfK[n-i] = k;
    }
    for(int i = 0; i < n; ++i)
        vectorLambdaNew[i] = matrixAk[i][i]; 
    return k;
}