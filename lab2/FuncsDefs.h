#pragma once
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#define my_type double
//#define my_type float	

void EprintEquation(my_type**& A, my_type*& B, const int& size, std::string s);
void EprintMatrix(my_type**& A, const int& size, std::string s);
void EprintVector(my_type*& B, const int& size, std::string s);

int EAllocateMemory(my_type**& A, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U,
    my_type*& B, my_type*& Bcopy, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& XR, my_type*& vectorD,
    my_type*& BF, my_type**& copyA, my_type**& DplusL, my_type**& DplusLI, my_type**& BI, my_type**& mBuffer,
    my_type**& T, my_type**& Q, my_type**& R, my_type*& copyB, my_type*& buffy, const int& n);
int EFreeMemory(my_type**& A, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U,
    my_type*& B, my_type*& Bcopy, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& XR, my_type*& vectorD,
    my_type*& BF, my_type**& copyA, my_type**& DplusL, my_type**& DplusLI, my_type**& BI, my_type**& mBuffer,
    my_type**& T, my_type**& Q, my_type**& R, my_type*& copyB, my_type*& buffy, const int& n);
int EReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type**& matrixA,
              my_type*& vectorB, my_type*& vectorBcopy, my_type**& copyA, my_type**& DplusL, my_type*& copyB, const int& n);

int EWriteVector(std::string fileNameOutput, const std::string label, my_type*& vector, const int& n);
int EWriteMatrix(const std::string fileNameOutput, const std::string label, my_type**& matrix, const int& n);
int EWriteData(std::string fileNameA, std::string fileNameB, std::string fileNameXZ, std::string fileNameXR, my_type**& matrixA,
              my_type*& vectorXZ, my_type*& vectorXR, my_type*& vectorB, const int& n);
int EWriteMatrixAB(const int& n);

my_type ECubicVectorNorm(my_type*& p, const int& size);
my_type EOctahedralVectorNorm(my_type*& p, const int& size);
//my_type ECubicMatrixNorm(my_type**& p, const int& size);
my_type EOctahedralMatrixNorm(my_type**& p, const int& size);
my_type EOctahedralMatrixNormClassic(my_type**& p, const int& size);

my_type*& EminusVectors(my_type*& X1, my_type*& X2, my_type*& BF, const int& size);

void ESwapVectors(my_type*& X1, my_type*& X2, const int& size);
void EMinusX(my_type*& X, const int& size);
int EZeidelMethod3D1(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& XZ2, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, const int& size);
my_type EResidue3d1(my_type**& matrixA, my_type*& vectorB, my_type*& vectorX, my_type*& vectorD, my_type*& buf, const int& size);

int EZeidelMethod3D2(my_type**& matrixA, my_type*& vectorB, my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, const int& size);
int ERelaxationMethod1(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& XR1, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, my_type& w, const int& size);
int ERelaxationMethod2(my_type**& matrixA, my_type*& vectorB, my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, my_type& w, const int& size);

my_type ECubicMatrixNorm(my_type**& p, const int& size);
void ECalculations(my_type**& matrixA, my_type*& vectorB, my_type*& vectorBcopy, my_type*& vectorX1, my_type*& vectorX2, my_type*& vectorXZ,
    my_type*& vectorXR, my_type*& matrixC_U, my_type*& matrixC_D, my_type*& matrixC_L, my_type*& vectorD, my_type*& BF, my_type& normC_L,
    my_type& normC_D, my_type& normC_U, my_type**& copyA, my_type**& DplusL, my_type**& DplusLI, my_type**& BI, my_type**& mBuffer, my_type**& T,
    my_type**& Q, my_type**& R, my_type*& copyB, my_type*& buffy, const int& size);
bool EMatrixDecomposition(my_type**& C, my_type*& matrixC_U, my_type*& matrixC_D, my_type*& matrixC_L,
                          my_type& normC_L, my_type& normC_D, my_type& normC_U, const int& n);
bool EisNotRight(my_type*& X1, my_type*& X2, my_type& norm, my_type*& BF, const int& size);
my_type ECountAccuracyZeidel(my_type**& matrixA, my_type*& vectorB, my_type**& DplusL, my_type**& DplusLI, my_type**& mBuffer, my_type**& T, my_type**& Q,
    my_type**& R, const int& size);
my_type ECountAccuracyRelaxation(my_type**& matrixA, my_type*& vectorB, my_type**& DplusL, my_type**& DplusLI, my_type**& mBuffer, my_type**& T, my_type**& Q,
    my_type**& R, my_type& w, const int& size);
my_type*& EplusVectors(my_type*& X1, my_type*& X2, my_type*& BF, const int& size);
my_type EMyCubicVectorNorm(my_type*& p, const int& size);

my_type EtypeResidue(my_type**& A, my_type*& B, my_type*& X, int norm_type, const int& size);
int EMatrixMult(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n);
void EMatrixTransform(my_type**& A, my_type*& B, my_type*& X, const int& size);
int EZeidelMethod3D(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& matrixC_L,
                   my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, const int& size);
int EZeidelMethod(my_type**& matrixC, my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type& norm,
                                                                            my_type*& BF, const int& size);
int ERelaxationMethod(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& matrixC_L,
                       my_type*& matrixC_D, my_type*& matrixC_U, my_type& norm, my_type*& BF, my_type& w,const int& size);

void EResidue3d(my_type**& matrixA, my_type*& vectorB, my_type*& vectorX, my_type*& BF, my_type*& buf, const int& size);

int EMatrixMult(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n);
//int EWriteMatrix(const std::string label, my_type**& matrix, const size_t n);
int EMatrixMultV(my_type**& matrixA, my_type**& matrixB, my_type**& matrixResult, const size_t n);
int EGetMatrixI(my_type**& matrix, const size_t n);
int EMatrixCopy(my_type**& matrixPaste, my_type**& matrixCopy, const size_t n);
int EMatrixTranspose(my_type**& matrixInit, my_type**& matrixResult, const size_t n);
int EMatrixInverse(my_type**& matrixA, my_type**& matrixInverted, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixBuffer, my_type*& vectorB, const size_t n);
int EMatrixInverseTR(my_type**& matrixT, my_type**& matrixR, my_type**& matrixInverted,
    my_type**& matrixBuffer, const size_t n);
int EQRDecomposer(my_type**& matrixA, my_type**& matrixT, my_type**& matrixQ,
    my_type**& matrixR, my_type*& vectorBuffer1, my_type*& vectorBuffer2, my_type*& vectorB, const size_t n);
int EQRCalculations(my_type**& matrixA, my_type**& matrixT, my_type**& matrixQ, my_type**& matrixR, my_type*& vectorB,
    my_type*& vectorX, my_type**& matrixBuffer1, my_type**& matrixBuffer2,
    my_type*& vectorBStarred, const size_t n);
int EReverseMotion(my_type**& matrixR, my_type*& vectorX, my_type*& vectorB, const size_t n);

void EZeidelMethod3DMarkCount(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XZ, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, int& iters, const int& size);
void ERelaxationMethodMarkCount(my_type*& B, my_type*& X1, my_type*& X2, my_type*& XR, my_type*& matrixC_L, my_type*& matrixC_D, my_type*& matrixC_U, my_type*& BF, my_type& w, int& iters, const int& size);






