#include "header.hpp"

void MatrixMult(double** &matrixLeft, double** &MatrixRight, double** &matrixResult, const size_t n)
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

void MatrixEDiff(double** &matrix, double** &matrixResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            if(i == j)
            {
                matrixResult[i] = 1 - matrix[i]; 
            }
            else
            {
                matrixResult[i] = - matrix[i];
            }
        }
    }
}









