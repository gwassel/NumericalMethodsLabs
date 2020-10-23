#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
#include "EmelinHeads.h"
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


    const std::string EpathData = "data/";

    std::string EfileNameMatrix;
    std::string EfileNameVector;
    std::string EfolderName;

    EfileNameMatrix = EpathData + "Ematrix1.txt";
    EfileNameVector = EpathData + "Evector1.txt";

    const std::string EfileMatrixEVecName = "res/matrixEVec1.txt";
    const std::string EfileMatrixEValName = "res/matrixEVal1.txt";
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