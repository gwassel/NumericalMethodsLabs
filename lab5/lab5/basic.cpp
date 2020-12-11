#include "Header.h"

void AllocateMemory(double*& mesh2D, double*& mesh3DX, double*& mesh3DY, double*& F, double**& u, double**& dF, double*& buf, double*& x0,
    double*& x1, double*& x2, int Fsize, int& meshSize2D, int& meshSize3DX, int& meshSize3DY) {

    mesh2D = new double[meshSize2D];
    mesh3DX = new double[meshSize3DX];
    mesh3DY = new double[meshSize3DY];

    F = new double[Fsize];

    u = new double* [Fsize];
    u[0] = new double[Fsize];
    u[1] = new double[Fsize];

    dF = new double* [Fsize];
    dF[0] = new double[Fsize];
    dF[1] = new double[Fsize];

    buf = new double[Fsize]{};
    x0 = new double[Fsize]{};
    x1 = new double[Fsize]{};
    x2 = new double[Fsize]{};
}



void CalculateOneEquation(double*& table, int& size, double A, double B, double (*f)(double)) {
    double epsilon = 1e-6;
    int t;
    std::vector<double> solsB;
    std::vector<double> solsN;

    //Bisection
    for (int i = 0; i < size - 1; i++) {
        t = HasSolution2D(table[i], table[i + 1], f);
        switch (t)
        {
        case 0:
            solsB.push_back(table[i]);
            //std::cout << "Found solution " << "[" << table[i] << "," << table[i + 1] << "]" << " " << solsB[solsB.size() - 1] << std::endl;
            break;
        case 1:
            solsB.push_back(BisectionMethod(table[i], table[i + 1], f));
            //std::cout << "Found solution " << "[" << table[i] << "," << table[i + 1] << "]" << " " << solsB[solsB.size() - 1] << std::endl;
            break;
        case 2:
            //std::cout << "No solution " << "[" << table[i] << "," << table[i + 1] << "]" << std::endl;
            break;
        }
    }
    //проверка последнего
    if (fabs(f(table[size - 1])) < epsilon) {
        solsB.push_back(table[size - 1]);
    }
    CheckForEqualElements(solsB);
    std::cout << "Solutions bisection(unique): " << std::endl;
    for (int i = 0; i < solsB.size(); i++) {
        std::cout << solsB[i] << std::endl;
    }

    //Newton method for one eq
    for (int i = 0; i < size - 1; i++) {
        t = HasSolution2D(table[i], table[i + 1], f);
        switch (t)
        {
        case 0:
            solsN.push_back(table[i]);
            //std::cout << "Found solution " << "[" << table[i] << "," << table[i + 1] << "]" << " " << solsN[solsN.size() - 1] << std::endl;
            break;
        case 1:
            solsN.push_back(NewtonOneEquation(table[i], table[i + 1], (table[i] + table[i + 1]) / 2.0, f));
            //std::cout << "Found solution " << "[" << table[i] << "," << table[i + 1] << "]" << " " << solsN[solsN.size() - 1] << std::endl;
            break;
        case 2:
            //std::cout << "No solution " << "[" << table[i] << "," << table[i + 1] << "]" << std::endl;
            break;
        }
    }
    if (fabs(f(table[size - 1])) < epsilon) {
        solsN.push_back(table[size - 1]);
    }
    CheckForEqualElements(solsN);
    std::cout << "Solutions newton(unique): " << std::endl;
    for (int i = 0; i < solsN.size(); i++) {
        std::cout << solsN[i] << std::endl;
    }
}