#include "header.hpp"

void ReadInit(const std::string path){}

void MAllocateMemory(){}

void AllocateMemory(double** &matrix, const size_t n)
{
    matrix = new double*[n];

    for(int i = 0; i < n; ++i)
    {
        matrix[i] = new double[n];
    }
}

void AllocateMemory(double* &vector, const size_t n)
{
    vector = new double[n];
}

void FreeMemory(){}

void FreeMemory(double** &matrix, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        delete[] matrix[i];
    }

    delete[] matrix;
}

void FreeMemory(double* &vector)
{
    delete[] vector;
}

void ReadData()
{
}

void ReadMatrix(std::string fileNameMatrix, double** &matrix, const size_t n)
{
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if ( !matrixFile.is_open() )
    {
        std::cerr << "Error: file with matrix is not open\n";
        return;
    }

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixFile >> matrix[i][j];
        }
    }

    matrixFile.close();
}

void ReadVector(std::string fileNameVector, double* &vector, const size_t n)
{
    std::ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if( !vectorFile.is_open() )
    {
        std::cerr << "Error: file with vector is not open\n"; 
        return;
    }

    for(int i = 0; i < n; ++i)
    {
        vectorFile >> vector[i];
    }
}

void WriteMatrixTeX(std::string fileOutputName, const double** &matrix, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileOutputName);

    fileOutput << "begin{matrix}\n";

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            fileOutput << matrix[i][j];
            if (j != n - 1)
                fileOutput << "& ";
        }
        fileOutput << "\\ \\ \n";
    }
}

