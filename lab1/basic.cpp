#include "funcs.hpp"

int ReadInit(const std::string path, std::string &fileNameA, std::string &fileNameB, std::string &fileNameQ,
        std::string &fileNameR, std::string &fileNameX, size_t &n)
{
    std::ifstream fileJsonInput;
    fileJsonInput.open(path);

    json objJson;
    fileJsonInput >> objJson;
    
    fileJsonInput.close();

    fileNameA = objJson["fileInputMatrixA"];
    fileNameB = objJson["fileInputVectorB"];
    fileNameQ = objJson["fileOutputMatrixQ"];
    fileNameR = objJson["fileOutputMatrixR"];
    fileNameX = objJson["fileOutputVectorX"];
    n = objJson["n"];

    return 0;
}

int AllocateMemory(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBuffer, const size_t n)
{
    AllocateMemory(matrixA, n);
    AllocateMemory(matrixT, n);
    AllocateMemory(matrixQ, n);
    AllocateMemory(matrixR, n);
    AllocateMemory(vectorB, n);
    AllocateMemory(vectorX, n);
    
    AllocateMemory(matrixBuffer1, n);
    AllocateMemory(matrixBuffer2, n);
    AllocateMemory(vectorBuffer, n);
    return 0;
}

int AllocateMemory(my_type** &matrix, const size_t n)
{
    matrix = new my_type*[n];

    for(int i = 0; i < n; ++i)
    {
        matrix[i] = new my_type[n];
    }
    return 0;
}

int AllocateMemory(my_type* &vector, const size_t n)
{
    vector = new my_type[n];

    return 0;
}

int ReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type** &matrixA,
        my_type* &vectorB, const size_t n)
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

int FreeMemory(my_type** &matrixA, my_type** &matrixT, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBuffer, const size_t n)
{
    FreeMemory(matrixA, n);
    FreeMemory(matrixT, n);
    FreeMemory(matrixQ, n);
    FreeMemory(matrixR, n);
    FreeMemory(vectorB, n);
    FreeMemory(vectorX, n);

    FreeMemory(matrixBuffer1, n);
    FreeMemory(matrixBuffer2, n);
    FreeMemory(vectorBuffer, n);

    return 0;
}

int FreeMemory(my_type** &matrix, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        delete[] matrix[i];
    }

    delete[] matrix;

    return 0;
}

int FreeMemory(my_type* &vector, const size_t n)
{
    delete[] vector;

    return 0;
}

int WriteMatrix(const std::string fileNameOutput, const std::string label, my_type** &matrix, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);
    
    fileOutput << label << "\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            fileOutput << matrix[i][j] << " ";
        }
        fileOutput << "\n";
    }
    fileOutput << "\n";
    
    fileOutput.close();

    return 0;
}

int WriteMatrix(const std::string label, my_type** &matrix, const size_t n)
{
    std::cout << label << "\n";

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            std::cout << matrix[i][j] << " ";    
        }
        std::cout << "\n";
    }
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

int WriteVector(const std::string label, my_type* &vector, const size_t n)
{
    std::cout << label << "\n";

    for(int i = 0; i < n; ++i)
    {
        std::cout << vector[i] << " ";
    }
    return 0;
}

int WriteJsonCfgsExample(const std::string path)
{
    json j = {
        {"n", 4},
        {"fileInputMatrixA", "data/matrixA"},
        {"fileInputVectorB", "data/vectorB"},
        {"fileOutputMatrixQ", "data/matrixQ"},
        {"fileOutputMatrixR", "data/matrixR"},
        {"fileOutputVectorX", "data/vectorX"}
    };

    std::ofstream fileJsonOutput;
    fileJsonOutput.open(path);

    fileJsonOutput << j << std::endl;
    fileJsonOutput.close();

    return 0;
}