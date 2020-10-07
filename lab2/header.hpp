#include <iostream>
#include <fstream>
#include <cmath>

#include "../nlohmann/json.hpp"
using json = nlohmann::json;

//basic functions
void ReadInit(const std::string, std::string&, std::string&, std::string&, std::string&, std::string&, std::string&, std::string&, std::string&, std::string&, std::string&, double &thau, double &epsilon, size_t &n);
void ReadData(const std::string, const std::string, const std::string, const std::string, double**&, double*&, double*&, const size_t);

void MAllocateMemory(double**&, double*&, double**&, double*&, double*&, double*&, const size_t);
void MAllocateBuffer(double*&, const size_t);
void MFreeMemory(double**&, double*&, double**&, double*&, double*&, double*&, const size_t);
void MFreeBuffer(double*&, const size_t);

void WriteMatrix(const std::string, double**&, const size_t);
void WriteVector(const std::string, double*&, const size_t);

void MWriteData(const std::string, double**&, double*&, double*&, const size_t);

void LogVector(std::string, double*&, const size_t);

//matrix functions
void MatrixMult(const double* const*, const double* const*, double**&, const size_t);
void MatrixMult(const double* const*, const double*, double*&, const size_t);
void MatrixCopy(double**&, const double* const*, const size_t);
void MatrixDot(double**&, const double, const size_t);

void VectorCopy(double*&, const double*, const size_t);
void VectorAdd(const double*, const double*, double*&, const size_t);
void VectorDiff(const double*, const double*, double*&, const size_t);
void VectorDot(double*&, const double, const size_t);

void MatrixEDiff(double**&, const size_t);

double CubicMatrixNorm(const double* const*, const size_t);
double CubicVectorNorm(const double*, const size_t);

//Molochkov functions
void MCalculations(const double* const*, const double*, double**&, double*&, double*&, double*&, double*&, const double, const double, const size_t);
