#include <iostream>
#include <fstream>
#include <string>
#include "../nlohmann/json.hpp"

#define my_type double

using json = nlohmann::json;

//header
int ReadInit(const std::string path); //read configs

int AllocateMemory(my_type** &matrixA, my_type** &matrixQ, my_type** matrixR, my_type* &vectorB, my_type* &columnX, const size_t n);

int ReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type** &matrixA, my_type* &vectorB, const size_t n); //read matrix and column

int Calculations(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, my_type* &columnX, const size_t n); 

int WriteData(std::string fileNameQ, std::string fileNameR, std::string fileNameX, 
        my_type** &matrixQ, my_type** &matrixR, my_type* &vectorX, const size_t n); //write results in output file


int FreeMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, my_type* &columnX, const size_t n); 

int MatrixMult(my_type** &matrixA,  my_type** &matrixB, const size_t n); //mult of 2 n*n matrix
int MatrixMunt(my_type** &matrix, my_type* &vector, const size_t n); //mult of n*n matrix on 1*n

int QRDecomposer(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, const size_t n);

int WriteMatrix(const std::string FileNameOutput, const std::string label, my_type** &matrix, const size_t n);
int WriteVector(const std::string FileNameOutput, const std::string label, my_type* &vector, const size_t n);

//init
int AllocateMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, const size_t n)
{
    matrixA = new my_type*[n];
    matrixQ = new my_type*[n];
    matrixR = new my_type*[n];

    for(int i = 0; i < n; ++i)
    {
        matrixA[i] = new my_type[n];
        matrixQ[i] = new my_type[n];
        matrixR[i] = new my_type[n];
    }
    
    vectorB = new my_type[n];
    
    return 0;
}

int ReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type** &matrixA, my_type* &vectorB, const size_t n)
{
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);
    
    if( !matrixFile.is_open() )
    {
        std::cerr << "Error: file with matrix is not open\n";
        return 1;
    }

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixFile >> matrixA[i][j];
        }
    }

    matrixFile.close();

    std::ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if( !vectorFile.is_open() )
    {
        std::cerr << "Error: file with vector is not open\n";
        return 1;
    }

    for(int i = 0; i < n; ++i)
    {
        vectorFile >> vectorB[i];
    }

    return 0;
}

int WriteData(std::string fileNameQ, std::string fileNameR, std::string fileNameX,
        my_type** &matrixQ, my_type** &matrixR, my_type* &vectorX, const size_t n)
{
    WriteMatrix(fileNameQ, "matrix Q", matrixQ, n);
    WriteMatrix(fileNameR, "matrix R", matrixR, n);
    WriteVector(fileNameX, "vector X", vectorX, n);

    return 0;
}

int FreeMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, const size_t n)
{
    delete[] vectorB;

    for(int i = 0; i < n; ++i)
    {
        delete[] matrixA[i];
        delete[] matrixQ[i];
        delete[] matrixR[i];
    }

    delete[] matrixA;
    delete[] matrixQ;
    delete[] matrixR;

    return 0;
}

int WriteMatrix(const std::string fileNameOutput, const std::string label, my_type** &matrix, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);
    
    fileOutput << label << "\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j){
            fileOutput << matrix[i][j] << " ";
        }
        fileOutput << "\n";
    }
    fileOutput << "\n";
    
    fileOutput.close();

    return 0;
}

int WriteVector(std::string fileNameOutput, const std::string label, my_type* &vector, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    fileOutput << label << "\n";
    for(int i = 0; i < n; ++i)
    {
        fileOutput << vector[i] << "\n";
    }
    fileOutput << "\n";

    fileOutput.close();

    return 0;
}

int main()
{
    my_type** matrixA;
    my_type* vectorB;
    my_type* vectorX;
    my_type** matrixQ;
    my_type** matrixR;
    
    const std::string pathConfig = "configs/config.json";
    const std::string pathData = "data/";
    
    std::string fileNameMatrix = pathData + "matrix";
    std::string fileNameVector = pathData + "vector";
    
    std::string fileNameQ = pathData + "outputQ";
    std::string fileNameR = pathData + "outputR";
    std::string fileNameX = pathData + "outputX";

    size_t n = 5;

    AllocateMemory(matrixA, matrixQ, matrixR, vectorB, n);

    ReadData(fileNameMatrix, fileNameVector, matrixA, vectorB, n);   
    
    WriteData(fileNameQ, fileNameR, fileNameX, matrixA, matrixA, vectorB, n);

    FreeMemory(matrixA, matrixQ, matrixR, vectorB, n);

    return 0;
}
