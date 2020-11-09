#include <iostream>
#include <cmath>
#include <fstream>
#include <string>


int ReadData(const std::string fileNameMatrix, const std::string fileNameEigenValsInit, double**& A, double*& lambda, const int size);
//�������� �������
int EAllocateMemory(double**& matrixA, double*& vectorX1, double*& vectorX2, double*& lambda, double*& vectorBstar,
    double*& vectorBuffer, double**& matrixBuffer1, double**& matrixBuffer2, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixEigenVectors, double**& matrixC,
    const int size);
void ECalculations(double**& matrixA, double*& vectorX1, double*& vectorX2, double*& lambda, double*& vectorBstar,
    double*& vectorBuffer, double**& matrixBuffer1, double**& matrixBuffer2, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixEigenVectors, double**& matrixC,
    const int size);
int EWriteData(std::string fileNameEVec, std::string fileNameEVal, double**& matrixEVec, double*& vectorEVals,
               const int& size);
int EFreeMemory(double**& matrixA, double*& vectorX1, double*& vectorX2, double*& lambda, double*& vectorBstar,
    double*& vectorBuffer, double**& matrixBuffer1, double**& matrixBuffer2, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixEigenVectors, double**& matrixC,
    const int size);

//Molochkov
void MAllocateMemory(double** &matrixH, double** &matrixAk, double** &matrixQ, double** &matrixR, double** &matrixBuffer, double* &vectorLambdaOld, 
        double* &vectorLambdaNew, const size_t n);
void MCalculations();
void MWriteData();
void MFreeMemory();
int SimpleQRIterations(double** &matrixAk, double** &matrixQ, double** &matrixR, double** &matrixBuffer, double* &vectorLambdaOld, 
        double* &vectorLambdaNew, double accuracy, const size_t n);
int ShiftQRIterations(double** &matrixAk, double** &matrixQ, double** &matrixR, double** &matrixBuffer, double* &vectorLambdaOld, 
        double* &vectorLambdaNew, double accuracy, const size_t n);

//������ � ����
int WriteVector(std::string fileNameOutput, double*& vector, const int& n);
int WriteMatrix(const std::string fileNameOutput,  double**& matrix, const int& n);

//�����
double CubicVectorNorm(double*& p, const int& size);
double CubicMatrixNorm(double**& p, const int& size);

//����� � �������
void PrintMatrix(double**& A, const int& size, std::string s);
void PrintVector(double*& B, const int& size, std::string s);

//��������� ������ � ��������
int MatrixMultVector(double**& matrix, double*& vector, double*& vectorResult, const size_t n);
int MatrixMultMatrix(double**& matrixA, double**& matrixB, double**& matrixResult, const size_t n);
double VectorMultVector(double*& vectorX1, double*& vectorX2, const size_t n);
void MatrixResE(double** &matrix, const size_t n, double c);

//QR
int MatrixMult(double**& matrix, double*& vector, double*& vectorResult, const size_t n);
int MatrixMultV(double**& matrixA, double**& matrixB, double**& matrixResult, const size_t n);
int GetMatrixI(double**& matrix, const size_t n);
int MatrixCopy(double**& matrixPaste, double**& matrixCopy, const size_t n);
int MatrixTranspose(double**& matrixInit, double**& matrixResult, const size_t n);
int MatrixInverse(double**& matrixA, double**& matrixInverted, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixBuffer, double*& vectorB, const size_t n);
int MatrixInverseTR(double**& matrixT, double**& matrixR, double**& matrixInverted,
    double**& matrixBuffer, const size_t n);
int QRDecomposer(double**& matrixA, double**& matrixT, double**& matrixQ,
    double**& matrixR, double*& vectorB, const size_t n);
int QRCalculations(double**& matrixA, double**& matrixT, double**& matrixQ, double**& matrixR, double*& vectorB,
    double*& vectorX, double**& matrixBuffer1, double**& matrixBuffer2,
    double*& vectorBStarred, const size_t n);
int ReverseMotion(double**& matrixR, double*& vectorX, double*& vectorB, const size_t n);
int ConditionNumberQR(double**& matrixR, double**& matrixT, double*& vector, const size_t column, const size_t n);
void QRDecomposerLite(double**& matrixA, double**& matrixT, double**& matrixQ, double**& matrixR, const size_t n);
int ConditionNumberQRLite(double**& matrixR, double**& matrixT, const size_t column, const size_t n);

//Hessenberg
void HessenbergForm(double** &matrixA, double** &matrixH, const size_t n);


void AllocateMatrix(double**&, const size_t);
void AllocateVector(double*&, const size_t);
void FreeMatrix(double*&, const size_t);
void FreeVector(double*&, const size_t);
