#include <iostream>
#include <fstream>
#include <cmath>

//basic functions
void ReadInit();
void ReadData(std::string, std::string, std::string, double**&, double*&, double*&, const size_t);

void MAllocateMemory(double**&, double*&, double**&, double*&, double*&, double*&, const size_t);
void MAllocateBuffer(double*&, const size_t);
void MFreeMemory(double**&, double*&, double**&, double*&, double*&, double*&, const size_t);
void MFreeBuffer(double*&, const size_t);

void WriteMatrix(const std::string, double**&, const size_t);
void WriteVector(const std::string, double*&, const size_t);

void MWriteData(std::string, std::string, std::string, double**&, double*&, double*&, const size_t);

void LogVector(std::string, double*&, const size_t);

//matrix functions
void MatrixMult(double**&, double**&, double**&, const size_t);
void MatrixMult(double**&, double*&, double*&, const size_t);
void MatrixCopy(double**&, double**&, const size_t);
void MatrixDot(double**&, const double, const size_t);

void VectorCopy(double*&, double*&, const size_t);
void VectorAdd(double*&, double*&, double*&, const size_t);
void VectorDiff(double*&, double*&, double*&, const size_t);
void VectorDot(double*&, const double, const size_t);

void MatrixEDiff(double**&, const size_t);

double CubicMatrixNorm(double**&, const size_t);
double CubicVectorNorm(double*&, const size_t);

//Molochkov functions
void MCalculations(double**&, double*&, double**&, double*&, double*&, double*&, double*&, const size_t);
