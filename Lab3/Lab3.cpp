#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
#include "EmelinHeads.h"
#include <nlohmann/json.hpp>
#define my_type double

int main()
{
    std::cout << "--------------------------------------------------------Emelin-------------------------------------------------------" << std::endl;
    const int Esize = 4;

    my_type** EmatrixA;
    my_type** EmatrixT;
    my_type** EmatrixQ;
    my_type** EmatrixR;
    my_type** EmatrixC;
    my_type** EmatrixBuffer1;
    my_type** EmatrixBuffer2;
    my_type** EmatrixEigenVectors;

    my_type* Elambda;
    my_type* EvectorX1;
    my_type* EvectorX2;
    my_type* EvectorBstar;
    my_type* EvectorBuffer;
    my_type* EvectorLambda;

    nlohmann::json Ej;
    std::ifstream Efile("menu/Econfig.json");
    const int ESYSTEM = 0;
    Efile >> Ej;

    const std::string EfileNameMatrix = Ej[ESYSTEM]["package_init"].get<std::string>() + Ej[ESYSTEM]["matrix_init_name"].get<std::string>();

    const std::string EfileMatrixEVecName = Ej[ESYSTEM]["package_res"].get<std::string>() + Ej[ESYSTEM]["eigen_vectors_file_res_name"].get<std::string>();
    const std::string EfileMatrixEValName = Ej[ESYSTEM]["package_res"].get<std::string>() + Ej[ESYSTEM]["eigen_values_file_res_name"].get<std::string>();
    
    //Убрать на случай другой системы лямбда
    EAllocateMemory(EmatrixA,EvectorX1,EvectorX2,EvectorLambda,EvectorBstar,EvectorBuffer,EmatrixBuffer1,
                    EmatrixBuffer2,EmatrixT,EmatrixQ,EmatrixR,EmatrixEigenVectors,EmatrixC,Esize);

    EReadData(EfileNameMatrix, EmatrixA, Esize);

    ECalculations(EmatrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);
    //Исправить, не тот массив 
    EWriteData(EfileMatrixEVecName, EfileMatrixEValName, EmatrixA, EvectorLambda, Esize);

    EFreeMemory(EmatrixA, EvectorX1, EvectorX2, EvectorLambda, EvectorBstar, EvectorBuffer, EmatrixBuffer1, EmatrixBuffer2,
        EmatrixT, EmatrixQ, EmatrixR, EmatrixEigenVectors, EmatrixC, Esize);

    system("pause");
}