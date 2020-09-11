#include <iostream>
#include <string>
#include "../nlohmann/json.hpp"

#define my_type double

using json = nlohmann::json;

int ReadInit(const std::string path); //read configs

int AllocateMemory(my_type** &matrixA, my_type* &columnB, const size_t n);

int ReadData(const std::string path, my_type** &matrixA, my_type* columnB, const size_t n); //read matrix and column

int Calculations(my_type** &matrixA, my_type* &matrixB);

int WriteData(my_type** matrixQ, my_type** matrixR, my_type* columnX); //write results in output file

int FreeMemory(my_type** &matrixA, my_type* &matrixB); 

int MatrixMult(my_type** &matrixA,  my_type** &matrixB, const size_t n); //mult of 2 n*n matrix
int MatrixMunt(my_type** &matrix, my_type* &column, const size_t n); //mult of n*n matrix on 1*n

int QRDecomposer(my_type** &matrixA, my_type** &matrixQ, my_type** &matrixR, const size_t n);


int AllocateMemory(my_type** &matrixA, my_type* &columnB, const size_t n)
{
    matrixA = new my_type*[n];
    
    for(int i = 0; i < n; ++i)
    {
        matrixA[i] = new my_type[n];
    }
    
    columnB = new my_type[n];
    
    return 0;
}

int ReadData(const std::string path, my_type** &matrixA, my_type* &columnB, const size_t n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            std::cin >> matrixA[i][j];
        }
    }

    for(int i = 0; i < n; ++i)
    {
        std::cin >> columnB[i];
    }

    return 0;
}


int main()
{
    my_type** matrixA;
    my_type* columnB;
    my_type* columnX;
    my_type** matrixQ;
    my_type** matrixR;
    
    const char pathInit[] = "configs/config.json";
    const char pathData[] = "data/";
    
    size_t n = 0;


    return 0;
}
