#include "Header.hpp"
const double epsilon = 1e-6;

int EAllocateMemory(double**& matrixA, double*& vectorX1, double*& vectorX2, double*& lambda, double*& vectorBstar,
                    double*& vectorBuffer, double**& matrixBuffer1, double**& matrixBuffer2, double**& matrixT,
                    double**& matrixQ, double**& matrixR, double**& matrixEigenVectors, double**& matrixC,
                    const int size) {
    matrixA = new double * [size];
    for (int i = 0; i < size; i++) {
        matrixA[i] = new double[size]{};
    }

    vectorX1 = new double[size]{};
    vectorX2 = new double[size]{};
    lambda = new double[size]{};

    vectorBstar = new double[size]{};

    vectorBuffer = new double[size]{};
    matrixBuffer1 = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixBuffer1[i] = new double[size]{};
    }
    matrixBuffer2 = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixBuffer2[i] = new double[size]{};
    }
    matrixT = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixT[i] = new double[size]{};
    }
    matrixQ = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixQ[i] = new double[size]{};
    }
    matrixR = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixR[i] = new double[size]{};
    }
    matrixEigenVectors = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixEigenVectors[i] = new double[size]{};
    }
    matrixC = new double * [size];
    for (int i = 0; i < size; ++i) {
        matrixC[i] = new double[size];
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrixC[i][j] = matrixA[i][j];
        }
    }

    return 0;
}
int ReadData(const std::string fileNameMatrix, const std::string fileNameEigenValsInit, double**& matrixA, double*& lambda, const int size) {
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with matrix is not open\n";
        return 1;
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrixFile >> matrixA[i][j];
        }
    }
    matrixFile.close();

    matrixFile.open(fileNameEigenValsInit);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with fileNameEigenValsInit is not open\n";
        return 1;
    }

    for (int i = 0; i < size; i++)
    {
            matrixFile >> lambda[i];
    }
    matrixFile.close();

    return 0;
}
void ECalculations(double**& matrixA, double*& vectorX1, double*& vectorX2, double*& lambda, double*& vectorBstar,
    double*& vectorBuffer, double**& matrixBuffer1, double**& matrixBuffer2, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixEigenVectors, double**& matrixC,
    const int size) {
    std::cout << "Epsilon: " << epsilon << std::endl;
    std::cout << std::endl;

    //����� �������� ����������� �����
    PrintVector(lambda, size, "Eigen values(init): ");
    double normY = 0.0;
    double lambdaPrev;
    int counter;
    //�������� 
    for (int k = 0; k < size; k++) {
        //������ �������� � ������������ ������
        //���� ������
        for (int i = 0; i < size; i++) {
            vectorX1[i] = 0.0;
        }
        vectorX1[k] = 1;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                matrixC[i][j] = matrixA[i][j];
            }
        }
        for (int i = 0; i < size; i++) {
            matrixC[i][i] -= lambda[k];
        }
        //������ �������
        QRCalculations(matrixC, matrixT, matrixQ, matrixR, vectorX1, vectorX2, matrixBuffer1, matrixBuffer2,
                        vectorBstar, size);

        //���������
        normY = 0.0;
        for (int i = 0; i < size; i++) {
            normY += vectorX2[i] * vectorX2[i];
        }
        normY = sqrt(normY);

        for (int i = 0; i < size; i++) {
            vectorX2[i] /= normY;
        }
        //������ ������� �������
        std::swap(vectorX1, vectorX2);
        counter = 1;
        //�� �� ����� � �����
        do{
            counter++;
            lambdaPrev = lambda[k];
            //���� ������������ �����
            MatrixMultVector(matrixA, vectorX1, vectorBuffer, size);
            lambda[k] = VectorMultVector(vectorBuffer, vectorX1, size);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    matrixC[i][j] = matrixA[i][j];
                }
            }
            for (int i = 0; i < size; i++) {
                matrixC[i][i] -= lambda[k];
            }

            QRCalculations(matrixC, matrixT, matrixQ, matrixR, vectorX1, vectorX2, matrixBuffer1,
                            matrixBuffer2, vectorBstar, size);

            normY = 0.0;
            for (int i = 0; i < size; i++) {
                normY += vectorX2[i] * vectorX2[i];
            }
            normY = sqrt(normY);

            for (int i = 0; i < size; i++) {
                vectorX2[i] /= normY;
            }
            std::swap(vectorX1, vectorX2);
        } while (fabs(lambdaPrev - lambda[k]) > epsilon);
        std::cout << "Iterations for lambda[" << k << "]: " << counter << std::endl;
        for (int i = 0; i < size; i++) {
            matrixEigenVectors[i][k] = vectorX1[i];
        }
    }
    std::cout << std::endl;
    PrintVector(lambda,size,"Eigen values(res): ");
    PrintMatrix(matrixEigenVectors,size,"Eigen vectors: ");
}
int EWriteData(std::string fileNameEVec, std::string fileNameEVal, double**& matrixEVec,
               double*& vectorEVals, const int& size) {

    WriteMatrix(fileNameEVec, matrixEVec, size);

    WriteVector(fileNameEVal, vectorEVals,size);

    return 0;
}
int EFreeMemory(double**& matrixA, double*& vectorX1, double*& vectorX2, double*& lambda, double*& vectorBstar,
    double*& vectorBuffer, double**& matrixBuffer1, double**& matrixBuffer2, double**& matrixT,
    double**& matrixQ, double**& matrixR, double**& matrixEigenVectors, double**& matrixC,
    const int size) {

    for (int i = 0; i < size; ++i)
    {
        delete[] matrixA[i];
    }
    delete[] matrixA;

    delete[] vectorX1;
    delete[] vectorX2;
    delete[] lambda;
    delete[] vectorBstar;
    delete[] vectorBuffer;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixBuffer1[i];
    }
    delete[] matrixBuffer1;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixBuffer2[i];
    }
    delete[] matrixBuffer2;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixT[i];
    }
    delete[] matrixT;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixQ[i];
    }
    delete[] matrixQ;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixR[i];
    }
    delete[] matrixR;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixEigenVectors[i];
    }
    delete[] matrixEigenVectors;
    for (int i = 0; i < size; ++i)
    {
        delete[] matrixC[i];
    }
    delete[] matrixC;
    return 0;
}
