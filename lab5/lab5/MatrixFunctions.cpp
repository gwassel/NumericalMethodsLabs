#include "Header.h"
//может выдавть минус ноль
void InverseMatrix(double**& matrix, double**& matrixInversed, int size) {
    double eps = 1e-6;
    double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    if (fabs(det) > eps) {
        matrixInversed[0][0] = matrix[1][1] / det;
        matrixInversed[0][1] = -matrix[0][1] / det;
        matrixInversed[1][0] = -matrix[1][0] / det;
        matrixInversed[1][1] = matrix[0][0] / det;
    }
    else {
        matrixInversed[0][0] = 0.0;
        matrixInversed[0][1] = 0.0;
        matrixInversed[1][0] = 0.0;
        matrixInversed[1][1] = 0.0;
        std::cout << "Det==0" << std::endl;
    }
}

void MatrixMultVector(const double* const* matrix, const double* vector, double*& vectorResult, const size_t n)
{
    double sum;
    for (int i = 0; i < n; ++i)
    {
        sum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }
}

double CubicMatrixNorm(const double* matrix, const size_t n)
{
    double sum = 0;
    double maxSum = 0;

    for (int i = 0; i < n; ++i)
    {
        maxSum += fabs(matrix[0 + i]);
        //std::cout << "Added " << i << std::endl;
    }

    for (int i = 1; i < n; ++i)
    {
        sum = 0;
        for (int j = 0; j < n; ++j)
        {
            sum += fabs(matrix[i * n + j]);
            //std::cout << "Added i: " << i<<", j: "<< j << std::endl;
        }
        if (sum > maxSum)
        {
            maxSum = sum;
        }
    }

    return sum;
}