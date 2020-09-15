#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../nlohmann/json.hpp"

#define my_type double

using json = nlohmann::json;

//header
int ReadInit(const std::string path, std::string &fileNameA, std::string &fileNameB, std::string &fileNameQ,
        std::string &fileNameR, std::string &fileNameX, size_t &n); //read configs

int AllocateMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBuffer, const size_t n);
int AllocateMemory(my_type** &matrix, const size_t n);
int AllocateMemory(my_type* &vector, const size_t n);

int ReadData(const std::string fileNameMatrix, const std::string fileNameVector, 
        my_type** &matrixA, my_type* &vectorB, const size_t n); //read matrix and column

int Calculations(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, 
        my_type* &vectorB, my_type* &columnX, my_type** &matrixBuffer1, 
        my_type** &matrixBuffer2, my_type* &vectorBStarred, const size_t n);

int WriteData(std::string fileNameQ, std::string fileNameR, std::string fileNameX, 
        my_type** &matrixQ, my_type** &matrixR, my_type* &vectorX, const size_t n); //write results in output file

int FreeMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, my_type* &columnX, const size_t n); 
int FreeMemory(my_type** &matrix, const size_t n);
int FreeMemory(my_type* &vector, const size_t n);

int MatrixMult(my_type** &matrixA,  my_type** &matrixB, my_type** &matrixResult, const size_t n); //mult of 2 n*n matrix
int MatrixMult(my_type** &matrix, my_type* &vector, my_type* &vectorResult, const size_t n); //mult of n*n matrix on 1*n
int MatrixCopy(my_type** &matrixPaste, my_type** &matrixCopy, const size_t n); // copy matrix
int MatrixTranspose(my_type** &matrixInit, my_type** &matrixResult, const size_t n); // transpose matrix init into matrix result
int VectorCopy(my_type* &vectorPaste, my_type* &vectorCopy, const size_t n);

int QRDecomposer(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type** &matrixT,
        my_type** &matrixBuffer, const size_t n);
int GetMatrixT(my_type** &matrixA, my_type** &matrixT, my_type** &matrixTi, const size_t n);
int GetMatrixI(my_type** &matrix, const size_t n);
int ReverseMotion(my_type** &matrixA, my_type* &vectorX, my_type* &vectorB, const size_t n); // solve Ax = b

int WriteMatrix(const std::string FileNameOutput, const std::string label, my_type** &matrix, const size_t n); // write matrix to file
int WriteMatrix(const std::string label, my_type** &matrix, const size_t n); // write matrix to console
int WriteVector(const std::string FileNameOutput, const std::string label, my_type* &vector, const size_t n);
int WriteVector(const std::string label, my_type* &vector, const size_t n);

int WriteJsonCfgsExample(const std::string path);

//init
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

int AllocateMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBuffer, const size_t n)
{
    AllocateMemory(matrixA, n);
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

int Calculations(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type* &vectorX, my_type** &matrixBuffer1, my_type** &matrixBuffer2, 
        my_type* &vectorBStarred, const size_t n)
{
    QRDecomposer(matrixA, matrixQ, matrixR, matrixBuffer1, matrixBuffer2, n);
    WriteMatrix("R", matrixR, n);
    MatrixMult(matrixR, vectorB, vectorBStarred, n);
    //b* = Tb    
    //Rx = b*
    ReverseMotion(matrixR, vectorX, vectorBStarred, n);

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

int FreeMemory(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type* &vectorB, 
        my_type** &matrixBuffer1, my_type** &matrixBuffer2, my_type* &vectorBuffer, const size_t n)
{
    FreeMemory(matrixA, n);
    FreeMemory(matrixQ, n);
    FreeMemory(matrixR, n);
    FreeMemory(vectorB, n);

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

int MatrixMult(my_type** &matrixA, my_type** &matrixB, my_type** &matrixResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            my_type sum = 0;
            for(int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}//matrixA * matrixB

int MatrixMult(my_type** &matrix, my_type* &vector, my_type* &vectorResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        my_type sum = 0;
        for(int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }

    return (0);
}//matrix * vector

int MatrixCopy(my_type** &matrixPaste, my_type** &matrixCopy, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixPaste[i][j] = matrixCopy[i][j];
        }
    }
    return 0;
} 

int MatrixTranspose(my_type** &matrixInit, my_type** &matrixResult, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixResult[i][j] = matrixInit[j][i];
        }
    }

    return 0;
}

int VectorCopy(my_type* &vectorPaste, my_type* &vectorCopy, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        vectorPaste[i] = vectorCopy[i];
    }

    return 0;
}

int QRDecomposer(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, my_type** &matrixT,
        my_type** &matrixBuffer, const size_t n)
{
    GetMatrixT(matrixA, matrixT, matrixBuffer, n);

    MatrixMult(matrixT, matrixA, matrixR, n); // R = TA
    // Q = transpose(T);
    MatrixTranspose(matrixT, matrixQ, n); // Q = transpose(T)

    return 0;
}

int GetMatrixT(my_type** &matrixA, my_type** &matrixT, my_type** &matrixTi, const size_t n)
{
    GetMatrixI(matrixT, n);
    GetMatrixI(matrixTi, n);

    my_type** matrixBuffer = new my_type*[n];

    for(int i = 0; i < n; ++i)
    {
        matrixBuffer[i] = new my_type[n];
    }
    

    
    for(int i = 0; i < n - 1; ++i)
    {
        for(int j = i + 1; j < n; ++j)
        {
            if(matrixA[i][j] != 0)
            {
                my_type c = matrixA[i][i];
                my_type s = matrixA[j][i];
                
                my_type radical = sqrt(c * c + s * s);
                
                c /= radical;
                s /= radical;
                
                matrixTi[i][i] = c;
                matrixTi[j][j] = c;
                matrixTi[j][i] = -s;
                matrixTi[i][j] = s;

                WriteMatrix("Ti", matrixTi, n);

                MatrixMult(matrixTi, matrixT, matrixBuffer, n);
                MatrixCopy(matrixT, matrixBuffer, n);
                std::cout << "i= " << i << "j= " << j << "c= " << c<< "s= " << s << "\n";
                WriteMatrix("T", matrixT, n);

                GetMatrixI(matrixTi, n);
            }
        }
    }



    for(int i = 0; i < n; ++i)
    {
        delete[] matrixBuffer[i];
    }

    delete[] matrixBuffer;

    return 0;
}

int GetMatrixR(my_type** &matrixR, my_type* &vectorBuffer1, my_type* &vectorBuffer2, const size_t n)
{

    for(int i = 0; i < n - 1; ++i)
    {
        for(int j = i + 1; j < n; ++j)
        {
            my_type c = matrixR[i][i];
            my_type s = matrixR[j][i];
                
            my_type radical = sqrt(c * c + s * s);
                
            c /= radical;
            s /= radical;
                
            for(int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixR[i][k] + s * matrixR[j][k]; //matrixR[i][k]
                vectorBuffer2[k] = (-s) * matrixR[i][k] + c * matrixR[j][k]; //matrixR[j][k]
            }

            for(int k = 0; k < n; ++k)
            {
                matrixR[i][k] = vectorBuffer1[k];
                matrixR[j][k] = vectorBuffer2[k];
            }
            WriteMatrix("R", matrixR, n);
        }
    }

    return 0;
}

int GetMatrixI(my_type** &matrix, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            if(i == j)
            {
                matrix[i][j] = 1;
            }
            else
            {
                matrix[i][j] = 0;
            }
        }
    }

    return 0;
}

int ReverseMotion(my_type** &matrixA, my_type* &vectorX, my_type* &vectorB, const size_t n)
{
    vectorX[0] = 1;
    for(int i = n - 1; i >= 0; --i)
    {
        my_type sum = 0;
        for(int j = i + 1; j < n; ++j)
        {
            sum += matrixA[i][j] * vectorX[j];
        }
        vectorX[i] = (vectorB[i] - sum) / matrixA[i][i];
    }

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

int main()
{
    my_type** matrixA = nullptr;
    my_type* vectorB = nullptr;
    my_type* vectorX = nullptr;
    my_type** matrixQ = nullptr;
    my_type** matrixR = nullptr;
    
    my_type** matrixBuffer1 = nullptr;
    my_type** matrixBuffer2 = nullptr;
    my_type* vectorBuffer = nullptr;

    const std::string pathConfig = "configs/config.json";
    const std::string pathData = "data/";
    
    //inputs
    std::string fileNameA = "";
    std::string fileNameB = "";
    
    //outputs
    std::string fileNameQ = "";
    std::string fileNameR = "";
    std::string fileNameX = "";
    
    size_t n = 0;

    ReadInit(pathConfig, fileNameA, fileNameB, fileNameQ, fileNameR, fileNameX, n);
    
    AllocateMemory(matrixA, matrixQ, matrixR, vectorB, vectorX, matrixBuffer1,
            matrixBuffer2, vectorBuffer, n);

    ReadData(fileNameA, fileNameB, matrixA, vectorB, n);   
   
    //Calculations(matrixA, matrixQ, matrixR, vectorB, vectorX, 
    //        matrixBuffer1, matrixBuffer2, vectorBuffer, n);

    // QRDecomposer(matrixA, matrixQ, matrixR, matrixBuffer1, matrixBuffer2, n);

    GetMatrixR(matrixA, vectorB, vectorX, n);

    WriteData(fileNameQ, fileNameR, fileNameX, matrixQ, matrixA, vectorX, n);

    FreeMemory(matrixA, matrixQ, matrixR, vectorB, matrixBuffer1,
            matrixBuffer2, vectorBuffer, n);

    WriteJsonCfgsExample(pathConfig);
    
    return 0;
}
