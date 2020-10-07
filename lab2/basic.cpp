#include "header.hpp"

void ReadInit(const std::string path, std::string &folderInput, std::string &fileNameMatrixA, std::string &fileNameVectorB, std::string &fileNameVectorX0, std::string &folderOutput, std::string &folderTest, std::string &fileTestJacobiA, std::string &fileTestJacobiB, std::string &fileTestTrueSolve, std::string &fileTestJacobiInitApprox, double &thau, double &epsilon, size_t &n)
{
    std::ifstream fileJsonInput;
    fileJsonInput.open(path);

    json objJson;
    fileJsonInput >> objJson;

    fileJsonInput.close();

    folderInput = objJson["folderInput"];
    fileNameMatrixA = objJson["fileNameMatrixA"];
    fileNameVectorB = objJson["fileNameVectorB"];
    fileNameVectorX0 = objJson["fileNameVectorX0"];

    folderOutput = objJson["folderOutput"];
    
    folderTest = objJson["folderTest"];
    fileTestJacobiA = objJson["fileTestJacobiMatrixA"];
    fileTestJacobiB = objJson["fileTestJacobiVectorB"];
    fileTestTrueSolve = objJson["fileTestJacobiTrueSolve"];
    fileTestJacobiInitApprox = objJson["fileTestJacobiInitApprox"];

    thau = objJson["thau"];
    epsilon = objJson["epsilon"];

    n = objJson["n"];
}

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

void MAllocateMemory(double** &matrixA, double* &vectorB, double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorY, const size_t n)
{
    AllocateMemory(matrixA, n);
    AllocateMemory(vectorB, n);
    AllocateMemory(matrixC, n);
    AllocateMemory(vectorXCurrent, n);
    AllocateMemory(vectorXFollow, n);
    AllocateMemory(vectorY, n);
}

void MAllocateBuffer(double* &vectorBuffer1, const size_t n)
{
    AllocateMemory(vectorBuffer1, n);
}

void FreeMemory(double** &matrix, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        delete[] matrix[i];
    }

    delete[] matrix;
}

void FreeMemory(double* &vector, const size_t n=0)
{
    delete[] vector;
}

void MFreeMemory(double** &matrixA, double* &vectorB, double** &matrixC, double* &vectorXCurrent, double* &vectorXFollow, double* &vectorY, const size_t n)
{
    FreeMemory(matrixA, n);
    FreeMemory(vectorB, n);
    FreeMemory(matrixC, n);
    FreeMemory(vectorXCurrent, n);
    FreeMemory(vectorXFollow, n);
    FreeMemory(vectorY, n);
}

void MFreeBuffer(double* &vectorBuffer1, const size_t n)
{
    FreeMemory(vectorBuffer1, n);
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

    int size = 0;

    matrixFile >> size;

    if( size != n)
    {
        std::cerr << "Error: matrix size is not correct\n";
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

    int size = 0;

    vectorFile >> size;

    if( size != n)
    {
        std::cerr << "Error: vector size is not correct\n";
    }

    for(int i = 0; i < n; ++i)
    {
        vectorFile >> vector[i];
    }
}

void ReadData(const std::string folder, const std::string fileNameMatrixA, const std::string fileNameVectorB, const std::string fileNameVectorX0, double** &matrixA, double* &vectorB, double* &vectorX0, const size_t n)
{
    std::string pathMatrixA = folder + fileNameMatrixA;
    std::string pathVectorB = folder + fileNameVectorB;
    std::string pathVectorX0 = folder + fileNameVectorX0;
    
    ReadMatrix(pathMatrixA, matrixA, n);
    ReadVector(pathVectorB, vectorB, n);
    ReadVector(pathVectorX0, vectorX0, n);
}

void WriteMatrix(const std::string fileNameOutput, double** &matrix, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
           fileOutput << matrix[i][j] << " "; 
        }
        fileOutput << "\n";
    }
}

void WriteVector(const std::string fileNameOutput, double* &vector, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for(int i = 0; i < n; ++i)
    {
        fileOutput << vector[i] << "\n";
    }
}

void MWriteData(const std::string folderOutput, double** &matrixC, double* &vectorX, double* &vectorY, const size_t n)
{
    std::string fileNameMatrixC = folderOutput + "matrixC";
    std::string fileNameVectorX = folderOutput + "vectorX";
    std::string fileNameVectorY = folderOutput + "vectorY";

    WriteMatrix(fileNameMatrixC, matrixC, n);
    WriteVector(fileNameVectorX, vectorX, n);
    WriteVector(fileNameVectorY, vectorY, n);
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

void LogVector(std::string label, double* &vector, const size_t n)
{
    std::cout << label << "\n";
    
    for(int i = 0; i < n; ++i)
    {
        std::cout << vector[i] << " ";
    }
    std::cout << "\n";
}
