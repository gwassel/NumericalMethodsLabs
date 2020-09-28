#include "header.hpp"

void MatrixMult(double** &matrixLeft, double** &matrixRight, double** &matrixResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            double sum = 0;
            for(int k = 0; k < n; ++k)
            {
                sum += matrixLeft[i][k] * matrixRight[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }
}

void MatrixMult(double** &matrix, double* &vector, double* &vectorResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        double sum = 0;
        for(int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }
}

void MatrixCopy(double** &matrixPaste, double** &matrixCopy, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixPaste[i][j] = matrixCopy[i][j];
        }
    }
}

void MatrixDot(double** &matrix, const double number, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrix[i][j] *= number;
        }
    }
}

void VectorDot(double* &vector, const double number, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        vector[i] *= number;
    }
}

void VectorCopy(double* &vectorPaste, double* &vectorCopy, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        vectorPaste[i] = vectorCopy[i];
    }
}

void VectorAdd(double* &vectorLeft, double* &vectorRight, double* &vectorResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        vectorResult[i] = vectorLeft[i] + vectorRight[i];
    }    
}

void VectorDiff(double* &vectorLeft, double* &vectorRight, double* &vectorResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        vectorResult[i] = vectorLeft[i] - vectorRight[i];
    } 
}

void MatrixEDiff(double** &matrix, const size_t n)
{
    std::cout << "matrixEDiff\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            if(i == j)
            {
                matrix[i][j] = 1 - matrix[i][j]; 
            }
            else
            {
                matrix[i][j] = - matrix[i][j];
            }
        }
    }
}

double CubicVectorNorm(double* &vector, const size_t n)
{
    double sum = 0;

    for(int i = 0; i < n; ++i)
    {
        sum += fabs(vector[i]);
    }

    return sum;
}

double CubicMatrixNorm(double** &matrix, const size_t n)
{
    double sum = 0;
    double maxSum = 0;

    for(int i = 0; i < n; ++i)
    {
        maxSum += fabs(matrix[0][i]);
    }
    
    for(int i = 1; i < n; ++i)
    {
        sum = 0;
        for(int j = 0; j < n; ++j)
        {
            sum += fabs(matrix[i][j]);
        }
        if(sum > maxSum)
        {
            maxSum = sum;
        }
    }
    
    return sum;
}







