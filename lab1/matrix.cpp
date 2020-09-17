#include "funcs.hpp"


int MatrixMult(my_type** &matrixA, my_type** &matrixB, my_type** &matrixResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            my_type sum = 0;
            for(int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}//matrixA * matrixB

int MatrixMult(my_type** &matrix, my_type* &vector, my_type* &vectorResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        my_type sum = 0;
        for(int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }

    return (0);
}//matrix * vector

int MatrixCopy(my_type** &matrixPaste, my_type** &matrixCopy, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixPaste[i][j] = matrixCopy[i][j];
        }
    }
    return 0;
} 

int MatrixTranspose(my_type** &matrixInit, my_type** &matrixResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixResult[i][j] = matrixInit[j][i];
        }
    }

    return 0;
}

int VectorCopy(my_type* &vectorPaste, my_type* &vectorCopy, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        vectorPaste[i] = vectorCopy[i];
    }

    return 0;
}

int GetMatrixT(my_type** &matrixA, my_type** &matrixT, my_type** &matrixTi, const size_t n)
{
    GetMatrixI(matrixT, n);
    GetMatrixI(matrixTi, n);

    my_type** matrixBuffer = new my_type*[n];

    for(int i = 0; i < n; ++i)
    {
        matrixBuffer[i] = new my_type[n];
    }
    
    for(int i = 0; i < n - 1; ++i)
    {
        for(int j = i + 1; j < n; ++j)
        {
            if(matrixA[i][j] != 0)
            {
                my_type c = matrixA[i][i];
                my_type s = matrixA[j][i];
                
                my_type radical = sqrt(c * c + s * s);
                
                c /= radical;
                s /= radical;
                
                matrixTi[i][i] = c;
                matrixTi[j][j] = c;
                matrixTi[j][i] = -s;
                matrixTi[i][j] = s;

                MatrixMult(matrixTi, matrixT, matrixBuffer, n);
                MatrixCopy(matrixT, matrixBuffer, n);

                GetMatrixI(matrixTi, n);
            }
        }
    }

    for(int i = 0; i < n; ++i)
    {
        delete[] matrixBuffer[i];
    }

    delete[] matrixBuffer;

    return 0;
}



int GetMatrixI(my_type** &matrix, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            if(i == j)
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