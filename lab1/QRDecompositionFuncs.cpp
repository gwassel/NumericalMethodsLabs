#include "funcs.hpp"

const double epsilon = 1e-12;

int ConditionNumber(my_type** &matrix, my_type* &vector, const size_t column, const size_t n)
{
    std::cout << "goes CONDNUM stage: " << column << std::endl;
    WriteMatrix("matrixR: ", matrix, n);
    size_t maxNumber = column;
    my_type maxValue = matrix[column][column];
    for(int i = column; i < n; ++i)
    {
        std::cout << "if " << fabs(matrix[i][column]) << " > " << fabs(maxNumber) << std::endl;
        if(fabs(matrix[i][column]) > fabs(maxValue))
        {
            maxValue = matrix[i][column];
            maxNumber = i;
        }
    }

    if(maxNumber != column) //if diagonal element is not max
    {
        std::cout << "stage " << column << " max number " << maxNumber << std::endl;
        std::swap(matrix[column], matrix[maxNumber]);
        std::swap(vector[column], vector[maxNumber]);
        WriteMatrix("matrixR: ", matrix, n);
    }
    else
    {
        std::cout << "maxValue:" << maxValue << std::endl;
        std::cout << "maxNumber:" << maxNumber << std::endl;
    }

    return 0;
}

int ReverseMotion(my_type** &matrixR, my_type* &vectorX, my_type* &vectorB, const size_t n)
{
    // if(abs(matrixR[n-1][n-1]) <= epsilon)
    // {
    //     std::cout << matrixR[n-1][n-1] << std::endl;
    //     std::cout << "Matrix is singular\n";
    //     return 1;
    // }
    // else{
        vectorX[0] = 1;
        for(int i = n - 1; i >= 0; --i)
        {
            my_type sum = 0;
            for(int j = i + 1; j < n; ++j)
            {
                sum += matrixR[i][j] * vectorX[j];
            }
            vectorX[i] = (vectorB[i] - sum) / matrixR[i][i];
        }
    // }
    return 0;     
}

int QRCalculations(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBStarred, const size_t n)
{
    MatrixCopy(matrixR, matrixA, n);
    //QRDecomposer2(matrixA, matrixQ, matrixR, matrixBuffer1, matrixBuffer2, n);
    QRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer1[0], matrixBuffer1[1], vectorB, n);
    MatrixMult(matrixT, vectorB, vectorBStarred, n);
    ReverseMotion(matrixR, vectorX, vectorBStarred, n);
    return 0;
}

int QRDecomposer(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, 
        my_type** &matrixR, my_type* &vectorBuffer1, my_type* &vectorBuffer2, my_type* &vectorB, const size_t n)
{
    MatrixCopy(matrixR, matrixA, n);

    GetMatrixI(matrixT, n);

    for(int i = 0; i < n - 1; ++i)
    {
        //тут выбор главного элемента
        ConditionNumber(matrixR, vectorB, i, n);

        for(int j = i + 1; j < n; ++j)
        {
            my_type c = matrixR[i][i];
            my_type s = matrixR[j][i];
                
            my_type radical = 1 / sqrt(c * c + s * s);
                
            c *= radical;
            s *= radical;
                
            for(int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixR[i][k] + s * matrixR[j][k]; //matrixR[i][k]
                vectorBuffer2[k] = (-s) * matrixR[i][k] + c * matrixR[j][k]; //matrixR[j][k]
            }

            for(int k = 0; k < n; ++k)
            {
                matrixR[i][k] = vectorBuffer1[k];
                matrixR[j][k] = vectorBuffer2[k];
            }

            for(int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixT[i][k] + s * matrixT[j][k];
                vectorBuffer2[k] = (-s) * matrixT[i][k] + c * matrixT[j][k];
            }

            for(int k = 0; k < n; ++k)
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

int MatrixInverseTR(my_type** &matrixT, my_type** &matrixR, my_type** &matrixInverted, 
        my_type** &matrixBuffer, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixBuffer[0][j] = matrixT[j][i];
        }

        ReverseMotion(matrixR, matrixInverted[i], matrixBuffer[0], n);
    }

    MatrixTranspose(matrixInverted, matrixBuffer, n);
    MatrixCopy(matrixInverted, matrixBuffer, n);

    return 0;
}

int MatrixInverse(my_type** &matrixA, my_type** &matrixInverted, my_type** &matrixT,
        my_type** &matrixQ, my_type** &matrixR, my_type** &matrixBuffer, my_type* &vectorB, const size_t n)
{
    QRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer[0], matrixBuffer[1], vectorB, n);

    MatrixInverseTR(matrixT, matrixR, matrixInverted, matrixBuffer, n);

    return 0;
}

