#include "Header.hpp"

void AllocateMatrix(double** &matrix, const size_t n)
{
    matrix = new double*[n];
    for(int i = 0; i < n; ++i)
        matrix[i] = new double[n];
}

void AllocateVector(double* &vector, const size_t n)
{
    vector = new double[n];
}

void FreeMatrix(double** &matrix, const size_t n)
{
    for(int i = 0; i < n; ++i)
        delete[] matrix[i];

    delete[] matrix;
}

void FreeVector(double* &vector, const size_t n)
{
    delete[] vector;
}


int WriteVector(std::string fileNameOutput, double*& vector, const int& n)
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

int WriteMatrix(const std::string fileNameOutput, double**& matrix, const int& n)
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

double CubicVectorNorm(double*& p, const int& size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(p[i]);
    }
    return sum;
}

double CubicMatrixNorm(double**& p, const int& size) {
    double sum;
    double maxSum = 0.0;
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

void PrintMatrix(double**& A, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void PrintVector(double*& B, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << B[i] << std::endl;
    }
    std::cout << std::endl;
}


int MatrixMultVector(double**& matrix, double*& vector, double*& vectorResult, const size_t n)
{
    double sum = 0.0;
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

int MatrixMultMatrix(double**& matrixA, double**& matrixB, double**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}

double VectorMultVector(double*& vectorX1, double*& vectorX2, const size_t size) {
    double res = 0.0;
    for (int i = 0; i < size; i++) {
        res += vectorX1[i] * vectorX2[i];
    }
    return res;
}


int MatrixMult(double**& matrix, double*& vector, double*& vectorResult, const size_t n)
{
    double sum = 0.0;
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

int MatrixMultV(double**& matrixA, double**& matrixB, double**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}

int GetMatrixI(double**& matrix, const size_t n)
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

int MatrixCopy(double**& matrixPaste, double**& matrixCopy, const size_t n)
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

int MatrixTranspose(double**& matrixInit, double**& matrixResult, const size_t n)
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

int MatrixInverse(double**& matrixA, double**& matrixInverted, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixBuffer, double*& vectorB, const size_t n)
{
    QRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer[0], matrixBuffer[1], vectorB, n);

    MatrixInverseTR(matrixT, matrixR, matrixInverted, matrixBuffer, n);

    return 0;
}

int MatrixInverseTR(double**& matrixT, double**& matrixR, double**& matrixInverted,
    double**& matrixBuffer, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixBuffer[0][j] = matrixT[j][i];
        }

        ReverseMotion(matrixR, matrixInverted[i], matrixBuffer[0], n);
    }

    MatrixTranspose(matrixInverted, matrixBuffer, n);
    MatrixCopy(matrixInverted, matrixBuffer, n);

    return 0;
}

int QRDecomposer(double**& matrixA, double**& matrixT, double**& matrixQ,
    double**& matrixR, double*& vectorBuffer1, double*& vectorBuffer2, double*& vectorB, const size_t n)
{
    MatrixCopy(matrixR, matrixA, n);

    GetMatrixI(matrixT, n);

    for (int i = 0; i < n - 1; ++i)
    {
        //��� ����� �������� ��������
        ConditionNumberQR(matrixR, matrixT, vectorB, i, n);

        for (int j = i + 1; j < n; ++j)
        {
            double c = matrixR[i][i];
            double s = matrixR[j][i];

            double radical = 1 / sqrt(c * c + s * s);

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

    MatrixTranspose(matrixT, matrixQ, n);

    return 0;
}

int QRCalculations(double**& matrixA, double**& matrixT, double**& matrixQ, double**& matrixR, double*& vectorB,
    double*& vectorX, double**& matrixBuffer1, double**& matrixBuffer2,
    double*& vectorBStarred, const size_t n)
{
    MatrixCopy(matrixR, matrixA, n);
    //QRDecomposer2(matrixA, matrixQ, matrixR, matrixBuffer1, matrixBuffer2, n);
    QRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer1[0], matrixBuffer1[1], vectorB, n);
    MatrixMult(matrixT, vectorB, vectorBStarred, n);
    ReverseMotion(matrixR, vectorX, vectorBStarred, n);
    return 0;
}

int ReverseMotion(double**& matrixR, double*& vectorX, double*& vectorB, const size_t n)
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
        double sum = 0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += matrixR[i][j] * vectorX[j];
        }
        vectorX[i] = (vectorB[i] - sum) / matrixR[i][i];
    }
    // }
    return 0;
}

int ConditionNumberQR(double**& matrixR, double**& matrixT, double*& vector, const size_t column, const size_t n)
{
    //std::cout << "goes CONDNUM stage: " << column << std::endl;
    //WriteMatrix("matrixR: ", matrixR, n);
    size_t maxNumber = column;
    double maxValue = matrixR[column][column];
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