#pragma once
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#define my_type double

//Основные функции
int EAllocateMemory(my_type**& matrixA, my_type*& vectorX1, my_type*& vectorX2, my_type*& lambda, my_type*& vectorBstar,
    my_type*& vectorBuffer, my_type**& matrixBuffer1, my_type**& matrixBuffer2, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixEigenVectors, my_type**& matrixC,
    const int size);
int EReadData(const std::string fileNameMatrix, my_type**& A, const int size);
void ECalculations(my_type**& matrixA, my_type*& vectorX1, my_type*& vectorX2, my_type*& lambda, my_type*& vectorBstar,
    my_type*& vectorBuffer, my_type**& matrixBuffer1, my_type**& matrixBuffer2, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixEigenVectors, my_type**& matrixC,
    const int size);
int EWriteData(std::string fileNameEVec, std::string fileNameEVal, my_type**& matrixEVec, my_type*& vectorEVals,
               const int& size);
int EFreeMemory(my_type**& matrixA, my_type*& vectorX1, my_type*& vectorX2, my_type*& lambda, my_type*& vectorBstar,
    my_type*& vectorBuffer, my_type**& matrixBuffer1, my_type**& matrixBuffer2, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixEigenVectors, my_type**& matrixC,
    const int size);

//Запись в файл
int EWriteVector(std::string fileNameOutput, my_type*& vector, const int& n);
int EWriteMatrix(const std::string fileNameOutput,  my_type**& matrix, const int& n);

//Нормы
my_type ECubicVectorNorm(my_type*& p, const int& size);
my_type ECubicMatrixNorm(my_type**& p, const int& size);

//Вывод в консоль
void EprintMatrix(my_type**& A, const int& size, std::string s);
void EprintVector(my_type*& B, const int& size, std::string s);

//Умножение матриц и векторов
int EMatrixMultVector(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n);
int EMatrixMultMatrix(my_type**& matrixA, my_type**& matrixB, my_type**& matrixResult, const size_t n);
my_type EVectorMultVector(my_type*& vectorX1, my_type*& vectorX2, const size_t n);

//Методы QR-разложения
int EMatrixMult(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n);
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
int EConditionNumberQR(my_type**& matrixR, my_type**& matrixT, my_type*& vector, const size_t column, const size_t n);






























