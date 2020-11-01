#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "EmelinHeads.h"
#define my_type double

int EAllocateMemory(my_type**& matrixA, my_type*& vectorX1, my_type*& vectorX2, my_type*& lambda, my_type*& vectorBstar,
                    my_type*& vectorBuffer, my_type**& matrixBuffer1, my_type**& matrixBuffer2, my_type**& matrixT,
                    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixEigenVectors, my_type**& matrixC,
                    const int size) {
    matrixA = new my_type * [size];
    for (int i = 0; i < size; i++) {
        matrixA[i] = new my_type[size]{};
    }

    vectorX1 = new my_type[size]{};
    vectorX2 = new my_type[size]{};
    lambda = new my_type[size]{};

    vectorBstar = new my_type[size]{};

    vectorBuffer = new my_type[size]{};
    matrixBuffer1 = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixBuffer1[i] = new my_type[size]{};
    }
    matrixBuffer2 = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixBuffer2[i] = new my_type[size]{};
    }
    matrixT = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixT[i] = new my_type[size]{};
    }
    matrixQ = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixQ[i] = new my_type[size]{};
    }
    matrixR = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixR[i] = new my_type[size]{};
    }
    matrixEigenVectors = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixEigenVectors[i] = new my_type[size]{};
    }
    matrixC = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        matrixC[i] = new my_type[size];
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrixC[i][j] = matrixA[i][j];
        }
    }

    return 0;
}
int EReadData(const std::string fileNameMatrix, const std::string fileNameEigenValsInit, my_type**& matrixA, my_type*& lambda, const int size) {
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with matrix is not open\n";
        return 1;
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrixFile >> matrixA[i][j];
        }
    }
    matrixFile.close();

    matrixFile.open(fileNameEigenValsInit);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with fileNameEigenValsInit is not open\n";
        return 1;
    }

    for (int i = 0; i < size; i++)
    {
            matrixFile >> lambda[i];
    }
    matrixFile.close();

    return 0;
}
void ECalculations(my_type**& matrixA, my_type*& vectorX1, my_type*& vectorX2, my_type*& lambda, my_type*& vectorBstar,
    my_type*& vectorBuffer, my_type**& matrixBuffer1, my_type**& matrixBuffer2, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixEigenVectors, my_type**& matrixC,
    const int size) {
    //Вывод исходных собственных чисел
    EprintVector(lambda, size, "Eigen values(init): ");
    my_type normY = 0.0;
    //Алгоритм 
    for (int k = 0; k < size; k++) {
        //Первая итерация с приближенным лямбда
        //Ищем лямбда
        for (int i = 0; i < size; i++) {
            vectorX1[i] = 0.0;
        }
        vectorX1[k] = 1;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                matrixC[i][j] = matrixA[i][j];
            }
        }
        for (int i = 0; i < size; i++) {
            matrixC[i][i] -= lambda[k];
        }
        //Решаем систему
        EQRCalculations(matrixC, matrixT, matrixQ, matrixR, vectorX1, vectorX2, matrixBuffer1, matrixBuffer2, vectorBstar, size);

        //Нормируем
        normY = 0.0;
        for (int i = 0; i < size; i++) {
            normY += vectorX2[i] * vectorX2[i];
        }
        normY = sqrt(normY);

        //Меняем местами вектора
        for (int i = 0; i < size; i++) {
            vectorX2[i] /= normY;
        }
        for (int i = 0; i < size; i++) {
            normY = vectorX2[i];
            vectorX2[i] = vectorX1[i];
            vectorX1[i] = normY;
        }

        //То же самое в цикле
        for (int i = 0; i < 100; i++) {
            //Ищем собстсвенное число
            EMatrixMultVector(matrixA, vectorX1, vectorBuffer, size);
            lambda[k] = EVectorMultVector(vectorBuffer, vectorX1, size);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    matrixC[i][j] = matrixA[i][j];
                }
            }
            for (int i = 0; i < size; i++) {
                matrixC[i][i] -= lambda[k];
            }

            EQRCalculations(matrixC, matrixT, matrixQ, matrixR, vectorX1, vectorX2, matrixBuffer1,
                            matrixBuffer2, vectorBstar, size);

            normY = 0.0;
            for (int i = 0; i < size; i++) {
                normY += vectorX2[i] * vectorX2[i];
            }
            normY = sqrt(normY);

            for (int i = 0; i < size; i++) {
                vectorX2[i] /= normY;
            }
            for (int i = 0; i < size; i++) {
                normY = vectorX2[i];
                vectorX2[i] = vectorX1[i];
                vectorX1[i] = normY;
            }
        }
        for (int i = 0; i < size; i++) {
            matrixEigenVectors[i][k] = vectorX1[i];
        }
    }
    EprintVector(lambda,size,"Eigen values(res): ");
    EprintMatrix(matrixEigenVectors,size,"Eigen vectors: ");
}
int EWriteData(std::string fileNameEVec, std::string fileNameEVal, my_type**& matrixEVec,
               my_type*& vectorEVals, const int& size) {

    EWriteMatrix(fileNameEVec, matrixEVec, size);

    EWriteVector(fileNameEVal, vectorEVals,size);

    return 0;
}
int EFreeMemory(my_type**& matrixA, my_type*& vectorX1, my_type*& vectorX2, my_type*& lambda, my_type*& vectorBstar,
    my_type*& vectorBuffer, my_type**& matrixBuffer1, my_type**& matrixBuffer2, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixEigenVectors, my_type**& matrixC,
    const int size) {

    for (int i = 0; i < size; ++i)
    {
        delete[] matrixA[i];
    }
    delete[] matrixA;

    delete[] vectorX1;
    delete[] vectorX2;
    delete[] lambda;
    delete[] vectorBstar;
    delete[] vectorBuffer;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixBuffer1[i];
    }
    delete[] matrixBuffer1;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixBuffer2[i];
    }
    delete[] matrixBuffer2;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixT[i];
    }
    delete[] matrixT;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixQ[i];
    }
    delete[] matrixQ;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixR[i];
    }
    delete[] matrixR;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixEigenVectors[i];
    }
    delete[] matrixEigenVectors;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixC[i];
    }
    delete[] matrixC;
    return 0;
}

int EWriteVector(std::string fileNameOutput, my_type*& vector, const int& n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for (int i = 0; i < n; ++i)
    {
        fileOutput << vector[i] << "\n";
    }
    fileOutput << "\n";

    fileOutput.close();

    return 0;
}

int EWriteMatrix(const std::string fileNameOutput, my_type**& matrix, const int& n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j) {
            fileOutput << matrix[i][j] << " ";
        }
        fileOutput << "\n";
    }
    fileOutput << "\n";

    fileOutput.close();

    return 0;
}

my_type ECubicVectorNorm(my_type*& p, const int& size) {
    my_type sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(p[i]);
    }
    return sum;
}

my_type ECubicMatrixNorm(my_type**& p, const int& size) {
    my_type sum;
    my_type maxSum = 0.0;
    for (int j = 0; j < size; j++) {
        maxSum += fabs(p[0][j]);
    }
    for (int i = 1; i < size; i++) {
        sum = 0.0;
        for (int j = 0; j < size; j++) {
            sum += fabs(p[i][j]);
        }
        if (sum > maxSum)
            maxSum = sum;
    }
    return maxSum;
}

void EprintMatrix(my_type**& A, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void EprintVector(my_type*& B, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << B[i] << std::endl;
    }
    std::cout << std::endl;
}


int EMatrixMultVector(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n)
{
    my_type sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        sum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }

    return 0;
}

int EMatrixMultMatrix(my_type**& matrixA, my_type**& matrixB, my_type**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            my_type sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}

my_type EVectorMultVector(my_type*& vectorX1, my_type*& vectorX2, const size_t size) {
    my_type res = 0.0;
    for (int i = 0; i < size; i++) {
        res += vectorX1[i] * vectorX2[i];
    }
    return res;
}


int EMatrixMult(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n)
{
    my_type sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        sum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }

    return 0;
}

int EMatrixMultV(my_type**& matrixA, my_type**& matrixB, my_type**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            my_type sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}

int EGetMatrixI(my_type**& matrix, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                matrix[i][j] = 1;
            }
            else
            {
                matrix[i][j] = 0;
            }
        }
    }

    return 0;
}

int EMatrixCopy(my_type**& matrixPaste, my_type**& matrixCopy, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixPaste[i][j] = matrixCopy[i][j];
        }
    }
    return 0;
}

int EMatrixTranspose(my_type**& matrixInit, my_type**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixResult[i][j] = matrixInit[j][i];
        }
    }

    return 0;
}

int EMatrixInverse(my_type**& matrixA, my_type**& matrixInverted, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixBuffer, my_type*& vectorB, const size_t n)
{
    EQRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer[0], matrixBuffer[1], vectorB, n);

    EMatrixInverseTR(matrixT, matrixR, matrixInverted, matrixBuffer, n);

    return 0;
}

int EMatrixInverseTR(my_type**& matrixT, my_type**& matrixR, my_type**& matrixInverted,
    my_type**& matrixBuffer, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixBuffer[0][j] = matrixT[j][i];
        }

        EReverseMotion(matrixR, matrixInverted[i], matrixBuffer[0], n);
    }

    EMatrixTranspose(matrixInverted, matrixBuffer, n);
    EMatrixCopy(matrixInverted, matrixBuffer, n);

    return 0;
}

int EQRDecomposer(my_type**& matrixA, my_type**& matrixT, my_type**& matrixQ,
    my_type**& matrixR, my_type*& vectorBuffer1, my_type*& vectorBuffer2, my_type*& vectorB, const size_t n)
{
    EMatrixCopy(matrixR, matrixA, n);

    EGetMatrixI(matrixT, n);

    for (int i = 0; i < n - 1; ++i)
    {
        //тут выбор главного элемента
        EConditionNumberQR(matrixR, matrixT, vectorB, i, n);

        for (int j = i + 1; j < n; ++j)
        {
            my_type c = matrixR[i][i];
            my_type s = matrixR[j][i];

            my_type radical = 1 / sqrt(c * c + s * s);

            c *= radical;
            s *= radical;

            for (int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixR[i][k] + s * matrixR[j][k]; //matrixR[i][k]
                vectorBuffer2[k] = (-s) * matrixR[i][k] + c * matrixR[j][k]; //matrixR[j][k]
            }

            for (int k = 0; k < n; ++k)
            {
                matrixR[i][k] = vectorBuffer1[k];
                matrixR[j][k] = vectorBuffer2[k];
            }

            for (int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixT[i][k] + s * matrixT[j][k];
                vectorBuffer2[k] = (-s) * matrixT[i][k] + c * matrixT[j][k];
            }

            for (int k = 0; k < n; ++k)
            {
                matrixT[i][k] = vectorBuffer1[k];
                matrixT[j][k] = vectorBuffer2[k];
            }
            //WriteMatrix("R", matrixR, n);
        }
    }

    EMatrixTranspose(matrixT, matrixQ, n);

    return 0;
}

int EQRCalculations(my_type**& matrixA, my_type**& matrixT, my_type**& matrixQ, my_type**& matrixR, my_type*& vectorB,
    my_type*& vectorX, my_type**& matrixBuffer1, my_type**& matrixBuffer2,
    my_type*& vectorBStarred, const size_t n)
{
    EMatrixCopy(matrixR, matrixA, n);
    //QRDecomposer2(matrixA, matrixQ, matrixR, matrixBuffer1, matrixBuffer2, n);
    EQRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer1[0], matrixBuffer1[1], vectorB, n);
    EMatrixMult(matrixT, vectorB, vectorBStarred, n);
    EReverseMotion(matrixR, vectorX, vectorBStarred, n);
    return 0;
}

int EReverseMotion(my_type**& matrixR, my_type*& vectorX, my_type*& vectorB, const size_t n)
{
    // if(abs(matrixR[n-1][n-1]) <= epsilon)
    // {
    //     std::cout << matrixR[n-1][n-1] << std::endl;
    //     std::cout << "Matrix is singular\n";
    //     return 1;
    // }
    // else{
    vectorX[0] = 1;
    for (int i = n - 1; i >= 0; --i)
    {
        my_type sum = 0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += matrixR[i][j] * vectorX[j];
        }
        vectorX[i] = (vectorB[i] - sum) / matrixR[i][i];
    }
    // }
    return 0;
}

int EConditionNumberQR(my_type**& matrixR, my_type**& matrixT, my_type*& vector, const size_t column, const size_t n)
{
    //std::cout << "goes CONDNUM stage: " << column << std::endl;
    //WriteMatrix("matrixR: ", matrixR, n);
    size_t maxNumber = column;
    my_type maxValue = matrixR[column][column];
    for (int i = column; i < n; ++i)
    {
        //std::cout << "if " << fabs(matrixR[i][column]) << " > " << fabs(maxNumber) << std::endl;
        if (fabs(matrixR[i][column]) > fabs(maxValue))
        {
            maxValue = matrixR[i][column];
            maxNumber = i;
        }
    }

    if (maxNumber != column) //if diagonal element is not max
    {
        //std::cout << "stage " << column << " max number " << maxNumber << std::endl;
        std::swap(matrixR[column], matrixR[maxNumber]);
        std::swap(matrixT[column], matrixT[maxNumber]);
        std::swap(vector[column], vector[maxNumber]);
    }
    else
    {
        //std::cout << "maxValue:" << maxValue << std::endl;
        //std::cout << "maxNumber:" << maxNumber << std::endl;
    }

    return 0;
}