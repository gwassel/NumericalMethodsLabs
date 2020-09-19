#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm> 

#include "../nlohmann/json.hpp"
using json = nlohmann::json;

#define my_type double

//basic funcs
int ReadInit(const std::string path, std::string &fileNameA, std::string &fileNameB, std::string &fileNameQ,
        std::string &fileNameR, std::string &fileNameX, std::string &fileNameAInv, std::string &fileNameA_Ainv, size_t &n); //read configs

int AllocateMemory(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, my_type** &matrixAInverted,
        my_type* &vectorBuffer, const size_t n);
int AllocateMemory(my_type** &matrix, const size_t n);
int AllocateMemory(my_type* &vector, const size_t n);

int ReadData(const std::string fileNameMatrix, const std::string fileNameVector, 
        my_type** &matrixA, my_type* &vectorB, const size_t n); //read matrix and column

int WriteData(std::string fileNameQ, std::string fileNameR, std::string fileNameX, std::string fileNameA_AInv, std::string fileNameAInv,
        my_type** &matrixQ, my_type** &matrixR, my_type* &vectorX, my_type** &matrixA_AInv, 
        my_type** &matrixAInv, const size_t n); //write results in output file

int FreeMemory(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, my_type** &matrixAInverted,
        my_type* &vectorBuffer, const size_t n); 
int FreeMemory(my_type** &matrix, const size_t n);
int FreeMemory(my_type* &vector, const size_t n);

int WriteMatrix(const std::string FileNameOutput, const std::string label, my_type** &matrix, const size_t n); // write matrix to file
int WriteMatrix(const std::string label, my_type** &matrix, const size_t n); // write matrix to console
int WriteVector(const std::string FileNameOutput, const std::string label, my_type* &vector, const size_t n);
int WriteVector(const std::string label, my_type* &vector, const size_t n);

int WriteJsonCfgsExample(const std::string path);

//matrix funcs

int MatrixMult(my_type** &matrixA,  my_type** &matrixB, my_type** &matrixResult, const size_t n); //mult of 2 n*n matrix
int MatrixMult(my_type** &matrix, my_type* &vector, my_type* &vectorResult, const size_t n); //mult of n*n matrix on 1*n
int MatrixCopy(my_type** &matrixPaste, my_type** &matrixCopy, const size_t n); // copy matrix
int MatrixTranspose(my_type** &matrixInit, my_type** &matrixResult, const size_t n); // transpose matrix init into matrix result
int MatrixInverseTR(my_type** &matrixT, my_type** &matrixR, my_type** &matrixInverted, 
        my_type** &matrixBuffer, const size_t n);
int VectorCopy(my_type* &vectorPaste, my_type* &vectorCopy, const size_t n);

//QRD funcs
int QRDecomposer(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, 
        my_type* &vectorBuffer1, my_type* &vectorBuffer2, my_type* &vectorB, const size_t n);
int QRDecomposer2(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type** &matrixT,
        my_type** &matrixBuffer, const size_t n);
int GetMatrixT(my_type** &matrixA, my_type** &matrixT, my_type** &matrixTi, const size_t n);
int GetMatrixI(my_type** &matrix, const size_t n);
int ReverseMotion(my_type** &matrixA, my_type* &vectorX, my_type* &vectorB, const size_t n); // solve Ax = b
int QRCalculations(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBStarred, const size_t n);

//gauss funcs
void GaussCalculations(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP, my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, const int& size);
int WriteParameters(std::string fileNameOutput, const my_type dis1, 
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx);
int WriteData(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP,my_type**& matrixA, my_type**& matrixIA,
    my_type*& vectorX1, my_type*& vectorX2, my_type*& vectorB, const my_type dis1,
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx, const int& n);
my_type markCondition(my_type*& X1, my_type*& X2, my_type*& B1, my_type*& B2, my_type& dx, my_type& db, const int& size);
my_type*& minusVectors(my_type*& X1, my_type*& X2, const int& size);
int FreeMemory(my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, const int& n);
std::tuple<my_type,my_type,my_type> discrepancy(my_type**& A, my_type*& B, my_type*& X, const int& size);
std::tuple<my_type, my_type> ConditionNumber(my_type**& A, my_type**& IA, const int& size);
my_type OctahedralMatrixNorm(my_type** p, const int& size);
my_type OctahedralVectorNorm(my_type* p, const int& size);
my_type CubicMatrixNorm(my_type** p, const int& size);
my_type CubicVectorNorm(my_type* p, const int& size);
bool isDegenerate(my_type**& A, const int& i, const int& size);
void SwapColomns(my_type**& A, const int& j1, const int& j2, const int& size);
void SwapRows(my_type**& A, my_type*& B, const int& i1, const int& i2);
std::tuple<int, int> SearchMax(my_type**& A, const int& currentMinor, const int& size);
void DiagonalizeEquation(my_type**& A, my_type*& B, my_type*& X, const int& size, std::vector<std::tuple<int, int>>& permutations);
bool MatrixIsPrepared(my_type**& A, my_type*& B, const int& i, std::vector<std::tuple<int, int>>& permutations, const int& size);
void GaussMethod(my_type** A1, my_type* B, my_type*& X, const int& size, bool& flag);
void printVector(my_type*& B, const int& size, std::string s);
void printMatrix(my_type**& A, const int& size, std::string s);
void printEquation(my_type**& A, my_type*& B, const int& size, std::string s);


my_type**& InverseMatrix(my_type** A1, const int& size);