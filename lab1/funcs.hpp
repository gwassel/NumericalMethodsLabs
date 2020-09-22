#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm> 

#include "../nlohmann/json.hpp"
using json = nlohmann::json;

#define my_type double

//basic funcs
int EMain(const int Esize, my_type** EA, my_type** EIA, my_type* EB, my_type* EBF, my_type* EdB, my_type* EX1, 
        my_type* EX2, my_type Econd, my_type Edx, my_type Edb, bool EnoProblem, std::tuple<my_type, my_type, my_type> Ed11,
        std::tuple<my_type, my_type> Ec11, const std::string EpathConfig, const std::string EpathData, std::string EfileNameMatrix,
        std::string EfileNameVector, const std::string EfileMatrixAName, const std::string EfileMatrixAIName, const std::string EfileVectorX1Name,
        const std::string EfileVectorX2Name, const std::string EfileVectorBName, const std::string EfileParamsName, std::string EfolderName);
int MMain(my_type** &matrixA, my_type* &vectorB, my_type* &vectorX, my_type** &matrixQ, my_type** &matrixR, 
        my_type** &matrixT, my_type** &matrixAInverted, my_type** &matrixBuffer1, my_type** &matrixBuffer2,
        my_type* &vectorBuffer, std::string fileNameA, std::string fileNameB, std::string fileNameQ, 
        std::string fileNameR, std::string fileNameX, std::string fileNameAInv, std::string fileNameA_AInv, const size_t n);

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
int GetMatrixT(my_type** &matrixA, my_type** &matrixT, my_type** &matrixTi, const size_t n);
int GetMatrixI(my_type** &matrix, const size_t n);
int ReverseMotion(my_type** &matrixA, my_type* &vectorX, my_type* &vectorB, const size_t n); // solve Ax = b
int QRCalculations(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBStarred, const size_t n);

//gauss funcs

int EAllocateMemory(my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, my_type*& BF, const int& n);
int EReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type**& matrixA, my_type*& vectorB, const int& n);

//Вывод матриц и уравнений
void EprintEquation(my_type**& A, my_type*& B, const int& size, std::string s);
void EprintMatrix(my_type**& A, const int& size, std::string s);
void EprintVector(my_type*& B, const int& size, std::string s);
int EWriteVector(std::string fileNameOutput, const std::string label, my_type*& vector, const int& n);
int EWriteParameters(std::string fileNameOutput, const my_type dis1,
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx);
int EWriteMatrix(const std::string fileNameOutput, const std::string label, my_type**& matrix, const int& n);
int EWriteData(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP, my_type**& matrixA, my_type**& matrixIA,
    my_type*& vectorX1, my_type*& vectorX2, my_type*& vectorB, const my_type dis1,
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx, const int& n);

//Метод Гаусса
void EGaussMethod(my_type**& A, my_type*& B, my_type*& X, const int& size, bool& flag);

//Подготовка матрицы к следующему шагу в методе Гаусса, т.е. установка главного элемента и проверка вырожденности
bool EMatrixIsPrepared(my_type**& A, my_type*& B, const int& i, std::vector<std::tuple<int, int>>& permutations, const int& size);

//Правильная расстановка элементов в векторе Х после всех перестановок столбцов
void EDiagonalizeEquation(my_type**& A, my_type*& B, my_type*& X, const int& size, std::vector<std::tuple<int, int>>& permutations);

//Поиск максимума
std::tuple<int, int> ESearchMax(my_type**& A, const int& currentMinor, const int& size);

//Перестановка строк/столбцов
void ESwapRows(my_type**& A, my_type*& B, const int& i1, const int& i2);
void ESwapColomns(my_type**& A, const int& j1, const int& j2, const int& size);

//Проверка вырожденности
bool EisDegenerate(my_type**& A, const int& i, const int& size);

//Нормы
my_type ECubicVectorNorm(my_type* &p, const int& size);
my_type EOctahedralVectorNorm(my_type* &p, const int& size);
my_type ECubicMatrixNorm(my_type** &p, const int& size);
my_type EOctahedralMatrixNorm(my_type** &p, const int& size);

std::tuple<my_type, my_type> EConditionNumber(my_type**& A, my_type**& IA, const int& size);

std::tuple<my_type, my_type, my_type> Ediscrepancy(my_type**& A, my_type*& B, my_type*& X, my_type*& BF, const int& size);

int EFreeMemory(my_type**& matrixA, my_type**& IA, my_type*& vectorB, my_type*& dB, my_type*& X1, my_type*& X2, my_type*& EBF, const int& n);

my_type*& EminusVectors(my_type*& X1, my_type*& X2, my_type*& BF, const int& size);

my_type EmarkCondition(my_type*& X1, my_type*& X2, my_type*& B1, my_type*& B2, my_type& dx, my_type& db, const int& size);

void ECalculations(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP, my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, my_type*& BF, const int& size);

int EMatrixMult(my_type**& matrixA, my_type**& matrixB, const int& n);

//my_type**& myInverse(my_type**& A1, const int& size);



int EMatrixMult(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n);
int EWriteMatrix(const std::string label, my_type**& matrix, const size_t n);
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
//int ConditionNumber(my_type**& matrix, my_type*& vector, const size_t column, const size_t n);