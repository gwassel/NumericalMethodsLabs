#pragma once
#include "FuncsDefs.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
#define my_type double
//#define my_type float
//Константы
const double epsilon = 1e-4;

void ECalculations(my_type**& matrixA, my_type*& vectorB, my_type*& vectorBcopy, my_type*& vectorX1, my_type*& vectorX2, my_type*& vectorXZ,
    my_type*& vectorXR, my_type*& matrixC_U, my_type*& matrixC_D, my_type*& matrixC_L, my_type*& vectorD, my_type*& BF, my_type& normC_L,
    my_type& normC_D, my_type& normC_U, my_type**& copyA, my_type**& DplusL, my_type**& DplusLI, my_type**& BufBuf, my_type**& mBuffer, my_type**& T,
    my_type**& Q, my_type**& R, my_type*& copyB, my_type*& buffy, const int& size) {

    my_type w = 1.0;
    my_type minW = w;

    int iterations;
    int minIterations;
    bool wantMark = true;

    EMatrixTransform(matrixA, vectorB, vectorX1, size);

    bool needMark = EMatrixDecomposition(matrixA, matrixC_L, matrixC_D, matrixC_U, normC_L, normC_D, normC_U, size);

    //Подсчет нормы через В1 и В2, но более экономично, только для Зейделя
    my_type** B1plusB2 = new my_type * [2];
    for (int i = 0; i < 2; i++) {
        B1plusB2[i] = new my_type[size]{};
    }

    B1plusB2[0][0] = matrixC_U[0];
    for (int i = 1; i < size - 1; i++) {
        B1plusB2[0][i] = matrixC_L[i];
        B1plusB2[1][i] = matrixC_U[i];
    }
    B1plusB2[0][size - 1] = matrixC_L[size - 2];
    my_type maxSum = 0.0;
    maxSum = fabs(B1plusB2[0][0]);
    for (int i = 1; i < size - 1; i++) {
        if (fabs(B1plusB2[0][i]) + fabs(B1plusB2[1][i]) > maxSum) {
            maxSum = fabs(B1plusB2[0][i]) + fabs(B1plusB2[1][i]);
        }
    }
    if (fabs(B1plusB2[1][size - 1]) > maxSum) {
        maxSum = fabs(B1plusB2[1][size - 1]);
    }
    my_type zeidelsnorm = epsilon * (1 - maxSum) / EMyCubicVectorNorm(matrixC_U, size);

    std::cout << "(1 - ||B1 + B2||) / ||B2|| * epsilon = " << zeidelsnorm << std::endl;

    my_type* tochnoe = new my_type[size]{};
    for (int i = 0; i < size; i++) {
        tochnoe[i] = 2 - (i % 2);
    }

    std::cout << "****************************Zeidel Method****************************" << std::endl;
    my_type normZeid = ECountAccuracyZeidel(copyA, vectorBcopy, DplusL, DplusLI, mBuffer, T, Q, R, size);

    int myMark;
    if (needMark || wantMark) {
        EMatrixMult(copyA, vectorX1, BF, size);
        std::cout << "CubicVectorNorm(minusVectors(vectorX1, BF, buffy, size), size) " << ECubicVectorNorm(EminusVectors(vectorX1, BF, buffy, size), size) << std::endl;
        myMark = (int)std::log(epsilon * (1 - ECubicMatrixNorm(DplusL, size)) / ECubicVectorNorm(EminusVectors(vectorX1, BF, buffy, size), size)) / std::log(ECubicMatrixNorm(DplusL, size));
        std::cout << "Mark Iterations Zeidel " << myMark << std::endl;
        for (int i = 0; i < size; ++i) {
            vectorX1[i] = 0.0;
            vectorX2[i] = 0.0;
        }
        EZeidelMethod3DMarkCount(vectorB, vectorX1, vectorX2, vectorXZ, matrixC_L, matrixC_D, matrixC_U, BF, myMark, size);
        EWriteVector("MarkCountZeidel.txt", "vector XZ", vectorXZ, size);
        for (int i = 0; i < size; i++) {
            buffy[i] = vectorXZ[i];
        }
    }
    my_type** cpA = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        cpA[i] = new my_type[size]{};
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cpA[i][j] = copyA[i][j];
        }
    }
    
    EZeidelMethod3D(vectorB, vectorX1, vectorX2, vectorXZ, matrixC_L, matrixC_D, matrixC_U, /*zeidelsnorm*/normZeid, BF, size);
    //EprintVector(vectorBcopy, size, "NINIJPOK");
    EResidue3d(cpA, vectorBcopy, vectorXZ, vectorD, BF, size);
    //1
    //my_type normZeidMark = ECubicVectorNorm(EminusVectors(tochnoe, vectorXZ, BF, size), size);
    //std::cout << "|VectorX^k - VectorX^(k+1)| " << normZeidMark << std::endl;
    
    //for (int i = 0; i < size; i++) {
    //    buffy[i] = 0.0;
    //    vectorX1[i] = 0.0;
    //    vectorX2[i] = 0.0;
    //}
    //EZeidelMethod3D1(vectorB, vectorX1, vectorX2,buffy, vectorXZ, matrixC_L, matrixC_D, matrixC_U, zeidelsnorm/*normZeid*/, BF, size);
    //normZeidMark = ECubicVectorNorm(EminusVectors(buffy, vectorXZ, BF, size), size);
    //std::cout << "|VectorX^k - VectorX*|1 " << normZeidMark << std::endl;

    //EZeidelMethod3D2(cpA, vectorBcopy, vectorB, vectorX1, vectorX2, buffy, matrixC_L, matrixC_D, matrixC_U, BF, size);
    //normZeidMark = ECubicVectorNorm(EminusVectors(buffy, vectorXZ, BF, size), size);
    //std::cout << "|AX-B|2 " << normZeidMark << std::endl;//1

    //2
    //EZeidelMethod(matrixA, vectorB, vectorX1, vectorX2, vectorXZ, zeidelsnorm/*normZeid*/, BF, size);
    //EtypeResidue(matrixA, vectorBcopy, vectorXR, 1, size);
    std::cout << "**************************Relaxation Method**************************" << std::endl;
    vectorX1[0] = vectorBcopy[0];
    for (int i = 1; i < size; i++) {
        vectorX1[i] = 0.0;
    }
    for (int i = 0; i < size; i++) {
        vectorX2[i] = 0.0;
    }
    my_type normRel = ECountAccuracyRelaxation(copyA, vectorB, DplusL, DplusLI, mBuffer, T, Q, R, w, size);
    //4
    /*if (needMark||wantMark) {
        EMatrixMult(copyA, vectorX1, BF, size);
        myMark = (int)std::log(epsilon * (1 - ECubicMatrixNorm(DplusL, size)) / ECubicVectorNorm(EminusVectors(vectorX1, BF, buffy, size), size)) / std::log(ECubicMatrixNorm(DplusL, size));
        std::cout << "Mark Iterations Relaxation = " << myMark << " for w = " << w << std::endl;
        ERelaxationMethodMarkCount(vectorB,vectorX1,vectorX2, vectorXR, matrixC_L,matrixC_D,matrixC_U, BF,  w, myMark,size);
        EWriteVector("MarkCountRelaxation.txt", "vector XR", vectorXZ, size);
        for (int i = 0; i < size; i++) {
            buffy[i] = vectorXR[i];
        }
    }*/
    
    minIterations = ERelaxationMethod(vectorB, vectorX1, vectorX2, vectorXR, matrixC_L, matrixC_D, matrixC_U, normRel, BF, w, size);
    
    //3
    //my_type normRelMark = ECubicVectorNorm(EminusVectors(tochnoe, vectorXR, BF, size), size);
    //std::cout <<"|VectorXR^k - VectorX^(k+1)| " << normRelMark << std::endl;

    ////EMinusX(vectorXR, size);
    ////std::cout << "|XXXXXXXxx2| " << ECubicVectorNorm(vectorXZ, size) << std::endl;
    ////std::cout << "|XXXXXXXxx12| " << ECubicVectorNorm(vectorXR, size) << std::endl;
    ////std::cout << "|XXXXXXXxx123| " <<ECubicVectorNorm(EminusVectors(vectorXZ, vectorXR, BF, size), size) << std::endl;
    //
    //// for (int i = 0; i < 190; i++) {
    ////    w = w + 0.01;
    ////    iterations = ERelaxationMethod(vectorB, vectorX1, vectorX2, vectorXR, matrixC_L, matrixC_D, matrixC_U, normRel, BF, w, size);
    ////    if (minIterations > iterations) {
    ////        minW = w;
    ////        minIterations = iterations;
    ////    }
    ////}

    //std::cout << "Min iterations in relaxation method: "<< minIterations <<" Min w: "<< minW << std::endl;
    //for (int i = 0; i < size; i++) {
    //    vectorD[i] = 0.0;
    //    BF[i] = 0.0;
    //}
    //
    //EResidue3d(cpA, vectorBcopy, vectorXR, vectorD, BF, size);

    //for (int i = 0; i < size; i++) {
    //    buffy[i] = 0.0;
    //    BF[i] = 0.0;
    //    vectorX1[i] = 0.0;
    //    vectorX2[i] = 0.0;
    //}

    //ERelaxationMethod1(vectorB, vectorX1, vectorX2,buffy, vectorXR, matrixC_L, matrixC_D, matrixC_U, normRel, BF, w, size);
    //my_type normRelaMark = ECubicVectorNorm(EminusVectors(buffy, vectorXR, BF, size), size);
    //std::cout << "|VectorXR^k - VectorX*|1 " << normRelaMark << std::endl;
    //
    //ERelaxationMethod2(cpA, vectorBcopy, vectorB, vectorX1, vectorX2, buffy, matrixC_L, matrixC_D, matrixC_U, BF, w, size);
    //normRelaMark = ECubicVectorNorm(EminusVectors(buffy, vectorXR, BF, size), size);
    //std::cout << "|AX-B|2 " << normRelaMark << std::endl;//2

    //EtypeResidue(matrixA, vectorBcopy, vectorXR, 1, size);
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

int EAllocateMemory(my_type**& A, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U,
                    my_type*& B, my_type*& Bcopy, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& XR, my_type*& vectorD,
                    my_type*& BF, my_type**& copyA, my_type**& DplusL, my_type**& DplusLI, my_type**& BI, my_type**& mBuffer,
                    my_type**& T, my_type**& Q, my_type**& R, my_type*& copyB, my_type*& buffy,  const int& n){
    A = new my_type * [n];
    matrixC_L = new my_type[n-1]{};
    matrixC_D = new my_type[n]{};
    matrixC_U = new my_type[n-1]{};
    for (int i = 0; i < n; ++i) {
        A[i] = new my_type[n]{};
    }
    copyB = new my_type[n]{};
    buffy = new my_type[n]{};
    copyA = new my_type * [n];
    for (int i = 0; i < n; i++) {
        copyA[i] = new my_type[n]{};
    }
    DplusL = new my_type * [n];
    for (int i = 0; i < n; i++) {
        DplusL[i] = new my_type[n]{};
    }
    DplusLI = new my_type * [n];
    for (int i = 0; i < n; i++) {
        DplusLI[i] = new my_type[n]{};
    }
    BI = new my_type * [n];
    for (int i = 0; i < n; i++) {
        BI[i] = new my_type[n]{};
    }
    mBuffer = new my_type * [n];
    for (int i = 0; i < n; i++) {
        mBuffer[i] = new my_type[n]{};
    }
    T = new my_type * [n];
    for (int i = 0; i < n; i++) {
        T[i] = new my_type[n]{};
    }
    Q = new my_type * [n];
    for (int i = 0; i < n; i++) {
        Q[i] = new my_type[n]{};
    }
    R = new my_type * [n];
    for (int i = 0; i < n; i++) {
        R[i] = new my_type[n]{};
    }

    BF = new my_type[n]{};
    vectorD = new my_type[n]{};
    B = new my_type[n]{};
    Bcopy = new my_type[n]{};
    X1 = new my_type[n]{};
    X2 = new my_type[n]{};
    XZ = new my_type[n]{};
    XR = new my_type[n]{};
    return 0;
}

void EprintEquation(my_type**& A, my_type*& B, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << B[i];
        std::cout << std::endl;
    }
    std::cout << std::endl;
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

int EFreeMemory(my_type**& A, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U,
    my_type*& B, my_type*& Bcopy, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& XR, my_type*& vectorD,
    my_type*& BF, my_type**& copyA, my_type**& DplusL, my_type**& DplusLI, my_type**& BI, my_type**& mBuffer,
    my_type**& T, my_type**& Q, my_type**& R, my_type*& copyB, my_type*& buffy, const int& n){

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;

    delete[] matrixC_L;
    delete[] matrixC_D;
    delete[] matrixC_U;
    delete[] B;
    delete[] Bcopy;
    delete[] copyB;
    delete[] X1;
    delete[] X2;
    delete[] XZ;
    delete[] XR;
    delete[] vectorD;
    delete[] BF;
    delete[] buffy;

    for (int i = 0; i < n; ++i)
    {
        delete[] copyA[i];
    }
    delete[] copyA;

    for (int i = 0; i < n; ++i)
    {
        delete[] DplusL[i];
    }
    delete[] DplusL;

    for (int i = 0; i < n; ++i)
    {
        delete[] DplusLI[i];
    }
    delete[] DplusLI;

    for (int i = 0; i < n; ++i)
    {
        delete[] BI[i];
    }
    delete[] BI;

    for (int i = 0; i < n; ++i)
    {
        delete[] mBuffer[i];
    }
    delete[] mBuffer;

    for (int i = 0; i < n; ++i)
    {
        delete[] T[i];
    }
    delete[] T;

    for (int i = 0; i < n; ++i)
    {
        delete[] Q[i];
    }
    delete[] Q;

    for (int i = 0; i < n; ++i)
    {
        delete[] R[i];
    }
    delete[] R;
    return 0;
}

int EReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type**& matrixA, my_type*& vectorB, 
    my_type*& vectorBcopy, my_type**& copyA, my_type**& DplusL, my_type*& copyB, const int& n)
{
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with matrix is not open\n";
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixFile >> matrixA[i][j];
        }
    }

    matrixFile.close();

    std::ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if (!vectorFile.is_open())
    {
        std::cerr << "Error: file with vector is not open\n";
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        vectorFile >> vectorB[i];
    }

    for (int i = 0; i < n; ++i)
    {
        vectorBcopy[i] = vectorB[i];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            copyA[i][j] = matrixA[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            DplusL[i][j] = matrixA[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        copyB[i] = vectorB[i];
    }

    return 0;
}

int EWriteData(std::string fileNameA, std::string fileNameB, std::string fileNameXZ, std::string fileNameXR,
    my_type**& matrixA, my_type*& vectorXZ, my_type*& vectorXR, my_type*& vectorB, const int& n)
{
    EWriteMatrix(fileNameA, "matrix A", matrixA, n);

    EWriteVector(fileNameB, "vector B", vectorB, n);

    EWriteVector(fileNameXZ+".txt", "vector XZ", vectorXZ, n);

    EWriteVector(fileNameXR + ".txt", "vector XR", vectorXR, n);
    return 0;
}

/*int EWriteMatrix(const std::string label, my_type**& matrix, const size_t n)
{
    std::cout << label << "\n";

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    return 0;
}*/

int EWriteMatrix(const std::string fileNameOutput, const std::string label, my_type**& matrix, const int& n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    fileOutput << label << "\n";
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

int EWriteVector(std::string fileNameOutput, const std::string label, my_type*& vector, const int& n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    //fileOutput << label << "\n";
    for (int i = 0; i < n; ++i)
    {
        fileOutput << vector[i] << "\n";
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

my_type EMyCubicVectorNorm(my_type*& p, const int& size) {
    my_type maxSum = 0.0;
    maxSum = fabs(p[0]);
    for (int i = 1; i < size - 1; i++) {
        if (fabs(p[i]) > maxSum) {
            maxSum = fabs(p[i]);
        }
    }
    return maxSum;
}

/*my_type ECubicMatrixNorm(my_type**& p, const int& size) {
    my_type sum;
    my_type maxSum = 0.0;
    for (int j = 0; j < 2; j++) {
        maxSum += fabs(p[0][j]);
    }
    for (int i = 1; i < size-1; i++) {
        sum = 0.0;
        for (int j = i-1; j <= i+1; j++) {
            sum += fabs(p[i][j]);
        }
        if (sum > maxSum)
            maxSum = sum;
    }
    sum = 0.0;
    for (int j = size - 2; j < size; j++) {
        sum += fabs(p[size-1][j]);
    }
    if (sum > maxSum)
        maxSum = sum;
    return maxSum;
}*/

my_type EOctahedralVectorNorm(my_type*& p, const int& size) {
    my_type max = fabs(p[0]);
    for (int i = 1; i < size; i++) {
        if (max < fabs(p[i]))
            max = fabs(p[i]);
    }
    return max;
}

my_type EOctahedralMatrixNorm(my_type**& p, const int& size) {
    my_type sum;
    my_type maxSum = 0.0;
    for (int i = 0; i <= 1; i++) {
        maxSum += fabs(p[i][0]);
    }
    for (int i = 1; i < size-1; i++) {
        sum = 0.0;
        for (int j = i-1; j <= i+1; j++) {
            sum += fabs(p[j][i]);
        }
        if (sum > maxSum)
            maxSum = sum;
    }
    sum = 0.0;
    for (int j = size - 2; j < size; j++) {
        sum += fabs(p[j][size-1]);
    }
    if (sum > maxSum)
        maxSum = sum;
    return maxSum;
}

my_type EOctahedralMatrixNormClassic(my_type**& p, const int& size) {
    my_type sum;
    my_type maxSum = 0.0;
    for (int i = 0; i < size; i++) {
        maxSum += fabs(p[i][0]);
    }
    for (int i = 1; i < size; i++) {
        sum = 0.0;
        for (int j = 0; j < size; j++) {
            sum += fabs(p[j][i]);
        }
        if (sum > maxSum)
            maxSum = sum;
    }
    return maxSum;
}

bool EMatrixDecomposition(my_type**& C, my_type*& matrixC_U, my_type*& matrixC_D, my_type*& matrixC_L,my_type& normC_L, my_type& normC_D, my_type& normC_U, const int& n) {
    matrixC_L[0] = C[1][0];
    matrixC_D[0] = C[0][0];
    
    for (int j = 1; j < n-1; j++) {
        matrixC_U[j-1] = C[j - 1][j];
        matrixC_D[j] = C[j][j]; 
        matrixC_L[j] = C[j + 1][j];
    }

    matrixC_U[n - 2] = C[n - 2][n - 1];
    matrixC_D[n - 1] = C[n - 1][n - 1];

    normC_L = EMyCubicVectorNorm(matrixC_L,n-1);
    normC_U = EMyCubicVectorNorm(matrixC_U,n-1);

   /* printMatrix(C,n,"C");
    printVector(matrixC_L, n-1, "matrixC_L");
    printVector(matrixC_D, n, "matrixC_D");
    printVector(matrixC_U, n-1, "matrixC_U");*/

    my_type sum = normC_U + normC_L;

    std::cout << "Cubic vector norm of C_L = " << normC_L << std::endl;
    std::cout << "Cubic vector norm of C_U = " << normC_U << std::endl;
    std::cout << "Sum = " << sum << std::endl;
    if (sum < 1) {
        return true;
    }else
        return false;
}

void EMatrixTransform(my_type**& A, my_type*& B, my_type*& X, const int& size) {
    my_type buf;

    for (int i = 0; i < size; i++) {
        buf = A[i][i];
        A[i][i] = 0.0;
        for (int j = 0; j < size; j++) {
            A[i][j] = (-1) * A[i][j] / buf;
        }
        B[i] = (-1) * B[i] / buf;
    }

    X[0] = B[0];
    for (int i = 1; i < size; i++) {
        X[i] = 0.0;
    }
}



my_type*& EminusVectors(my_type*& X1, my_type*& X2, my_type*& BF, const int& size) {
    for (int i = 0; i < size; i++) {
        BF[i] = X1[i] - X2[i];
    }
    return BF;
}

my_type*& EplusVectors(my_type*& X1, my_type*& X2, my_type*& BF, const int& size) {
    for (int i = 0; i < size; i++) {
        BF[i] = X1[i] + X2[i];
    }
    return BF;
}

bool EisNotRight(my_type*& X1, my_type*& X2,my_type& norm, my_type*& BF, const int& size) {
    //std::cout << " Mistake norm: " << ECubicVectorNorm(EminusVectors(X1, X2, BF, size), size)<<std::endl;
    return (ECubicVectorNorm(EminusVectors(X1, X2,BF, size), size) > norm);
}

int EZeidelMethod3D(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U,my_type& norm, my_type*& BF, const int& size) {
    int k = 0;
    //std::ofstream fileOutput;
    //fileOutput.open("iterationsZeidel3D.txt"/*, std::ios_base::app*/);
    do {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        k++;
        X2[0] = B[0] + matrixC_U[0] * X1[1];
        /*for (int j = 0; j < size; ++j)
        {
            fileOutput << -X2[j] << " ";
        }
        fileOutput << "\n";*/
        for (int i = 1; i <= size-2; i++){
            X2[i] = B[i] + matrixC_L[i - 1] * X2[i - 1] + matrixC_U[i] * X1[i + 1];

            /*for (int j = 0; j < size; ++j)
            {
                fileOutput << -X2[j] << " ";
            }
            fileOutput << "\n";*/

        }
        X2[size - 1] = B[size - 1] + matrixC_L[size - 2] * X2[size - 2];

        /*for (int i = 0; i < size; ++i)
        {
            fileOutput << -X2[i] << " ";
        }
        fileOutput << "\n";*/

        //ESwapVectors(X1, X2, size);

    } while (EisNotRight(X1,X2,norm,BF, size));
    //while (ECubicVectorNorm(EminusVectors(X1, X2, BF, size), size) > epsilon);
    //fileOutput.close();
    std::cout << "Iterations in Zeidel method(3 diagonals): " << k << std::endl;
    EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XZ[i] = X2[i];
    }
    return k;
}

int EZeidelMethod3D1(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& XZ2, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U,my_type& norm, my_type*& BF, const int& size) {
    int k = 0;
    EMinusX(XZ2, size);
    do {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        k++;
        X2[0] = B[0] + matrixC_U[0] * X1[1];
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = B[i] + matrixC_L[i - 1] * X2[i - 1] + matrixC_U[i] * X1[i + 1];
        }
        X2[size - 1] = B[size - 1] + matrixC_L[size - 2] * X2[size - 2];
        std::cout << ECubicVectorNorm(EminusVectors(X2, XZ2, BF, size), size) << std::endl;
    } //while (EisNotRight(X1,X2,norm,BF, size));
    while (ECubicVectorNorm(EminusVectors(X2, XZ2, BF, size), size) > epsilon);
    //fileOutput.close();
    std::cout << "ITERATIONS1: " << k << std::endl;
    EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XZ[i] = X2[i];
    }
    EMinusX(XZ2, size);
    return k;
}


int EZeidelMethod3D2(my_type**& matrixA, my_type*& vectorB, my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, const int& size) {
    int k = 0;
    double res;
    for (int i = 0; i < size; i++) {
        X1[i] = 0.0;
        X2[i] = X1[i];
    }
    my_type* matrixD = new my_type[size]{};
    my_type* buf1 = new my_type[size]{};

    do {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        k++;
        X2[0] = B[0] + matrixC_U[0] * X1[1];
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = B[i] + matrixC_L[i - 1] * X2[i - 1] + matrixC_U[i] * X1[i + 1];
        }
        X2[size - 1] = B[size - 1] + matrixC_L[size - 2] * X2[size - 2];
        EMinusX(X2, size);
        res = EResidue3d1(matrixA,vectorB,X2,matrixD,buf1,size);
        EMinusX(X2, size);
    }// while (EisNotRight(X1,X2,norm,BF, size));
    while (res > epsilon);
    //fileOutput.close();
    std::cout << "ITERATIONSZ2: " << k << std::endl;
    EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XZ[i] = X2[i];
    }
    return k;
}

void EZeidelMethod3DMarkCount (my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, int& iters, const int& size) {

    for (int d = 0; d<iters; d++){

        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        } 
        X2[0] = B[0] + matrixC_U[0] * X1[1];
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = B[i] + matrixC_L[i - 1] * X2[i - 1] + matrixC_U[i] * X1[i + 1];

        }
        X2[size - 1] = B[size - 1] + matrixC_L[size - 2] * X2[size - 2];

    }
    EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XZ[i] = X2[i];
    }
}


int EZeidelMethod(my_type**& matrixC, my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type& norm, my_type*& BF, const int& size) {
    int k = 0;
    my_type temp;
    std::ofstream fileOutput;
    fileOutput.open("iterationsZeidel1.txt"/*, std::ios_base::app*/);
    do {
        for (int i = 0; i < size; i++) {
            X1[i] =  X2[i];
        }
        k++;
        for (int i = 0; i < size; i++) {
            temp = 0.0;

            for (int j = 0; j < i; j++) {
                temp += matrixC[i][j] * X2[j];
            }

            for (int j = i + 1; j < size; j++) {
                temp += matrixC[i][j] * X1[j];
            }

            X2[i] = B[i] + temp;

            for (int i = 0; i < size; ++i)
            {
                fileOutput << -X2[i] << " ";
            }
            fileOutput << "\n";
        }
        //ESwapVectors(X1, X2, size);


        //for (int i = 0; i < size; ++i)
        //{
        //    fileOutput << -X2[i] << " ";
        //}
        //fileOutput << "\n";

    } //while (EisNotRight(X1, X2, norm, BF, size));
    while (ECubicVectorNorm(EminusVectors(X1,X2,BF,size),size) > epsilon);
    fileOutput.close();
    std::cout << "Iterations in Zeidel method(matrix): " << k << std::endl;
    EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XZ[i] = X2[i];
    }
    return k;
}

my_type ECountAccuracyZeidel(my_type**& matrixA, my_type*& vectorB, my_type**& DplusL, my_type**& DplusLI, my_type**& mBuffer, my_type**& T,my_type**& Q,
    my_type**& R, const int& size) {
    EWriteMatrix("MatrDplusL.txt", "DplusL", DplusL, size);
    EWriteMatrix("MatrmatrixA.txt", "matrixA", matrixA, size);
    EMatrixInverse(DplusL, DplusLI, T, Q, R, mBuffer, vectorB, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j > i) {
                mBuffer[i][j] = -matrixA[i][j];
            }else
                mBuffer[i][j] = 0.0;
        }
    }

    EMatrixMultV(DplusLI, mBuffer,  DplusL, size);

    /*for (int i = 0; i < size; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                DplusL[i][j] = 1 - DplusL[i][j];
            }
            else {
                DplusL[i][j] = -DplusL[i][j];
            }
        }
    }*/

    EWriteMatrix("MATRIX_CZ.txt", "matrix C", DplusL, size);

    my_type normB = EOctahedralMatrixNormClassic(DplusL, size);
    std::cout << "Octahedral matrix norm (D + L) = " << normB << std::endl;
    normB = ((1 - normB) / normB) * epsilon;
    std::cout << "Octahedral accuracy (D + L) = " << normB << std::endl;

    normB = ECubicMatrixNorm(DplusL, size);
    std::cout << "Cubic matrix norm (D + L) = " << normB << std::endl;
    normB = ((1 - normB) / normB) * epsilon;
    std::cout << "Cubic accuracy (D + L) = " << normB << " (We take this)"<< std::endl;
    return normB;
}


//void ERelaxationMethod(my_type*& B, my_type*& X1, my_type*& X2, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, const int& size) {
//    int k = 0;
//    my_type w = -1.05;
//    do {
//        k++;
//        X2[0] = (1 - w) * X1[0] - w * (matrixC_U[0] * X1[1]) + w * B[0];
//
//        for (int i = 1; i <= size - 2; i++) {
//            X2[i] = (-w) * (matrixC_L[i - 1] * X2[i - 1]) + (1 - w) * X1[i] - w * (matrixC_U[i] * X1[i+1]) + w * B[i];
//        }
//        X2[size - 1] = (-w) * (matrixC_L[size - 2] * X2[size - 2]) + (1 - w) * X1[size - 1] +  w * B[0];
//        
//        SwapVectors(X1, X2, size);
//    } while (EisNotRight(X1, X2, norm, BF, size));
//    std::cout << "Its: " << k << std::endl;
//    MinusX(X1, size);
//}

int ERelaxationMethod(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, my_type& w, const int& size) {
    int k = 0;
    //std::ofstream fileOutput;
    //fileOutput.open("iterationsRelax.txt"/*, std::ios_base::app*/);
    //EprintVector(matrixC_D,size-1,"BBBBBBB");
    do {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        k++;
        X2[0] = (1 - w) * X1[0] + w * (matrixC_U[0] * X1[1]) - w * B[0];
        //for (int i = 0; i < size; ++i)
        //{
        //    fileOutput << X2[i] << " ";
        //}
        //fileOutput << "\n";
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = w * (matrixC_L[i - 1] * X2[i - 1]) + (1 - w) * X1[i] + w * (matrixC_U[i] * X1[i + 1]) - w * B[i];
            //for (int j = 0; j < size; ++j)
            //{
            //    fileOutput << X2[j] << " ";
            //}
            //fileOutput << "\n";
        }
        X2[size - 1] = w * (matrixC_L[size - 2] * X2[size - 2]) + (1 - w) * X1[size - 1] - w * B[size - 1];
        //for (int j = 0; j < size; ++j)
        //{
        //    fileOutput << X2[j] << " ";
        //}
        //fileOutput << "\n";
        //ESwapVectors(X1, X2, size);
    } //while (EisNotRight(X1, X2, norm, BF, size));
    while (ECubicVectorNorm(EminusVectors(X1, X2, BF, size), size) > epsilon);
    std::cout << "Iterations(Relaxation): " << k << std::endl;
    //fileOutput.close();
    //EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XR[i] = X2[i];
    }
    return k;
}

int ERelaxationMethod1(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& XR1, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, my_type& w, const int& size) {
    int k = 0;
    do {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        k++;
        X2[0] = (1 - w) * X1[0] + w * (matrixC_U[0] * X1[1]) - w * B[0];
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = w * (matrixC_L[i - 1] * X2[i - 1]) + (1 - w) * X1[i] + w * (matrixC_U[i] * X1[i + 1]) - w * B[i];
        }
        X2[size - 1] = w * (matrixC_L[size - 2] * X2[size - 2]) + (1 - w) * X1[size - 1] - w * B[size - 1];
    } //while (EisNotRight(X1, X2, norm, BF, size));
    while (ECubicVectorNorm(EminusVectors(X2, XR1, BF, size), size) > epsilon);
    std::cout << "Iterations(Relaxation1): " << k << std::endl;
    for (int i = 0; i < size; i++) {
        XR[i] = X2[i];
    }
    return k;
}


int ERelaxationMethod2(my_type**& matrixA, my_type*& vectorB, my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, my_type& w, const int& size) {
    int k = 0;
    double res;
    for (int i = 0; i < size; i++) {
        X1[i] = 0.0;
        X2[i] = X1[i];
    }

    my_type* matrixD = new my_type[size]{};
    my_type* buf1 = new my_type[size]{};
    do {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        k++;
        X2[0] = (1 - w) * X1[0] + w * (matrixC_U[0] * X1[1]) - w * B[0];
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = w * (matrixC_L[i - 1] * X2[i - 1]) + (1 - w) * X1[i] + w * (matrixC_U[i] * X1[i + 1]) - w * B[i];
        }
        X2[size - 1] = w * (matrixC_L[size - 2] * X2[size - 2]) + (1 - w) * X1[size - 1] - w * B[size - 1];
        //EMinusX(X2, size);
        res = EResidue3d1(matrixA, vectorB, X2, matrixD, buf1, size);
        //EMinusX(X2, size);
    } while (res > epsilon);
    std::cout << "Iterations(Relaxation): " << k << std::endl;
    //EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XR[i] = X2[i];
    }
    return k;
}


void ERelaxationMethodMarkCount(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, my_type& w,int& iters, const int& size) {
    for (int d = 0; d < iters;d++) {
        for (int i = 0; i < size; i++) {
            X1[i] = X2[i];
        }
        X2[0] = (1 - w) * X1[0] + w * (matrixC_U[0] * X1[1]) + w * B[0];
        for (int i = 1; i <= size - 2; i++) {
            X2[i] = w * (matrixC_L[i - 1] * X2[i - 1]) + (1 - w) * X1[i] + w * (matrixC_U[i] * X1[i + 1]) + w * B[i];
        }
        X2[size - 1] = w * (matrixC_L[size - 2] * X2[size - 2]) + (1 - w) * X1[size - 1] + w * B[0];
    }
    EMinusX(X2, size);
    for (int i = 0; i < size; i++) {
        XR[i] = X2[i];
    }
}

my_type ECountAccuracyRelaxation(my_type**& matrixA, my_type*& vectorB, my_type**& DplusL, my_type**& DplusLI, my_type**& mBuffer, my_type**& T, my_type**& Q,
    my_type**& R, my_type& w, const int& size) {

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j < i) {
                DplusL[i][j] = w * matrixA[i][j];
            }
            else if (i == j) {
                DplusL[i][j] = matrixA[i][j];
            }
            else {
                DplusL[i][j] = 0.0;
            }
        }
    }
    EMatrixInverse(DplusL, DplusLI, T, Q, R, mBuffer, vectorB, size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j < i) {
                DplusL[i][j] = w * matrixA[i][j];
            }
            else if (i == j) {
                DplusL[i][j] = matrixA[i][j];
            }
            else {
                DplusL[i][j] = 0.0;
            }
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            mBuffer[i][j] = DplusL[i][j] - w * matrixA[i][j];
        }
    }

    EMatrixMultV(DplusLI, mBuffer, DplusL, size);
    /*for (int i = 0; i < size; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                DplusL[i][j] = 1 - DplusL[i][j];
            }
            else {
                DplusL[i][j] = -DplusL[i][j];
            }
        }
    }*/


    my_type normB = EOctahedralMatrixNormClassic(DplusL, size);
    std::cout << "Octahedral matrix norm (D + w*L) = " << normB << std::endl;
    normB = ((1 - normB) / normB) * epsilon;
    std::cout << "Octahedral accuracy (D + w*L) = " << normB << std::endl;

    normB = ECubicMatrixNorm(DplusL, size);
    std::cout << "Cubic matrix norm (D + w*L) = " << normB << std::endl;
    normB = ((1 - normB) / normB) * epsilon;
    std::cout << "Cubic accuracy (D + w*L) = " << normB << " (We take this) " << std::endl;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (DplusL[i][j]<epsilon)
              DplusL[i][j] = 0.0;
        }
    }

    EWriteMatrix("MATRIX_CR.txt", "matrix C", DplusL, size);
    return normB;
}


void ESwapVectors(my_type*& X1, my_type*& X2,const int& size) {
    std::swap(X1,X2);
}

void EMinusX(my_type*& X, const int& size) {
    for (int i = 0; i < size; i++) {
        X[i] *= (-1);
    }
}

int EWriteMatrixAB(const int& n)
{
    std::ofstream fileOutput;
    //fileOutput.open("data/Ematrix1");
    fileOutput.open("data/Ematrix22.txt");
    /*for (int i = 0; i < n; i++)
    {
        //for (int j = 0; j < n; j++) {
        //    if (j == i) {
        //        fileOutput << j << " ";
        //    } else if (j == i-1||j == i+1) {
        //        fileOutput << 1 << " ";
        //    } else {
        //        fileOutput << 0 << " ";
        //    }


            if (j == i) {
                fileOutput << 8 << " ";
            }
            else if (j == i - 1) {
                fileOutput << 1 << " ";
            }
            else if (j == i + 1) {
                fileOutput << 6 << " ";
            }
            else {
                fileOutput << 0 << " ";
            }
        }
        fileOutput << "\n";
    }*/
    fileOutput << 5 << " ";
    fileOutput << 1 << " ";
    fileOutput << "\n";

    fileOutput << 3 << " ";
    fileOutput << 2 << " ";
    fileOutput << "\n";

    fileOutput.close();

    //fileOutput.open("data/Evector1");
    fileOutput.open("data/Evector22.txt");
    fileOutput << 6 << " ";
    fileOutput << "\n";
    fileOutput << 1 << " ";
    fileOutput << "\n";
    /*//fileOutput << 6 << "\n";
    fileOutput << 1 << "\n";
    for (int i = 1; i < n-1; i++)
    {
        //fileOutput << 10-2*(i%2) << "\n";
        fileOutput << i+1 << "\n";
    }
    //fileOutput << 9-3*(n%2) << "\n";
    fileOutput << n-1 << "\n";*/
    fileOutput.close();


    return 0;
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

my_type EtypeResidue(my_type**& A, my_type*& B, my_type*& X,int norm_type, const int& size) {

    my_type* buf = new my_type[size];
    my_type* buf2 = new my_type[size];

    EMatrixMult(A, X, buf, size);
    my_type norm;
    if (norm_type == 1)
        norm = ECubicVectorNorm(EminusVectors(B, buf, buf2, size),size);
    if (norm_type == 2)
        norm = EOctahedralVectorNorm(EminusVectors(B, buf, buf2, size),size);
    return norm;
}

void EResidue3d(my_type**& matrixA, my_type*& vectorB, my_type*& vectorX, my_type*& vectorD, my_type*& buf, const int& size) {
    vectorD[0] = matrixA[0][0] * vectorX[0] + matrixA[0][1] * vectorX[1];
    for (int i = 1; i <= size - 2; i++) {
        vectorD[i] = matrixA[i][i - 1] * vectorX[i-1] + matrixA[i][i] * vectorX[i] + matrixA[i][i+1] * vectorX[i+1];
    } 
    vectorD[size-1] = matrixA[size - 1][size - 2] * vectorX[size - 2] + matrixA[size - 1][size - 1] * vectorX[size - 1];    
    my_type normD = ECubicVectorNorm(EminusVectors(vectorD,vectorB,buf,size),size); 
    std::cout <<"Norm residue (Cubic) = "<< normD << std::endl;
    //normD = OctahedralVectorNorm(minusVectors(vectorD, vectorB, buf, size), size);
    //std::cout << "Norm D (Octahedral) = " << normD << std::endl;
}

my_type EResidue3d1(my_type**& matrixA, my_type*& vectorB, my_type*& vectorX, my_type*& vectorD, my_type*& buf, const int& size) {
    vectorD[0] = matrixA[0][0] * vectorX[0] + matrixA[0][1] * vectorX[1];
    for (int i = 1; i <= size - 2; i++) {
        vectorD[i] = matrixA[i][i - 1] * vectorX[i - 1] + matrixA[i][i] * vectorX[i] + matrixA[i][i + 1] * vectorX[i + 1];
    }
    vectorD[size - 1] = matrixA[size - 1][size - 2] * vectorX[size - 2] + matrixA[size - 1][size - 1] * vectorX[size - 1];
    my_type normD = ECubicVectorNorm(EminusVectors(vectorD, vectorB, buf, size), size);
    std::cout << "Norm residue (Cubic) = " << normD << std::endl;
    //normD = OctahedralVectorNorm(minusVectors(vectorD, vectorB, buf, size), size);
    //std::cout << "Norm D (Octahedral) = " << normD << std::endl;
    return normD;
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

int EMatrixInverse(my_type**& matrixA, my_type**& matrixInverted, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixBuffer, my_type*& vectorB, const size_t n)
{
    EQRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer[0], matrixBuffer[1], vectorB, n);

    EMatrixInverseTR(matrixT, matrixR, matrixInverted, matrixBuffer, n);

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
}//matrixA * matrixB



