#include "Header.h"

double ChordMethodIteration(double A, double x, double (*f)(double)) {
    return x - f(x) * (A - x) / (f(A) - f(x));
}

double NewtonOneEquation(double A, double B, double x, double (*f)(double)) {
    double x1 = x;
    int iters = 1;
    double x2 = x1 - f(x1) / Derivative(f, x1);
    double eps = 1e-6;
    int count = 0;
    const int N = 100;

    while (fabs(x1 - x2) > eps) {
        x2 = x1 - f(x1) / Derivative(f, x1);
        if (x2 > A || x2 < B) {
            if (Derivative(f, (A + B) / 2) * SecondDerivative(f, (A + B) / 2) >= 0) {
                count++;
                for (int i = 0; i < count; i++) {
                    x2 = ChordMethodIteration(B, x2, f);
                }
            }
            else if (Derivative(f, (A + B) / 2) * SecondDerivative(f, (A + B) / 2) < 0) {
                count++;
                for (int i = 0; i < count; i++) {
                    x2 = ChordMethodIteration(A, x2, f);
                }
            }
            if (count > N) {
                break;
                std::cout << "BAD" << std::endl;
            }
        }
        std::swap(x1, x2);
        iters++;
    }
    std::cout <<"Newt "<< iters << std::endl;
    return x2;
}

//если внутри есть разрыв то ошибется
double BisectionMethod(double A, double B, double (*f)(double)) {
    int iters = 1;
    double a = A;//>0
    double b = B;//<0
    double epsilon = 1e-6;
    double c = (a + b) / 2;
    double fc = f(c);
    //&& fabs(b-a) > epsilon???
    while (fabs(fc) > epsilon && fabs(b - a) > epsilon) {
        c = (a + b) / 2;
        if ((f(a) * f(c)) > 0) a = c;
        else b = c;
        fc = f(c);
        iters++;
    }
    std::cout <<"Bi "<< iters << std::endl;
    return c;
}


int NewtonSystem(double*& F, double**& dF, double**& u, double*& buf, double*& x0, double*& x1, double*& x2,
    double (*f1)(double, double), double (*f2)(double, double), int size) {
    double eps = 1e-6;
    double norm;
    int iters = 0;

    x1[0] = x0[0];
    x1[1] = x0[1];

    F[0] = f1(x1[0], x1[1]);
    F[1] = f2(x1[0], x1[1]);

    dF[0][0] = PartialDerivativeX(f1, x1[0], x1[1]);
    dF[0][1] = PartialDerivativeY(f1, x1[0], x1[1]);
    dF[1][0] = PartialDerivativeX(f2, x1[0], x1[1]);
    dF[1][1] = PartialDerivativeY(f2, x1[0], x1[1]);

    //RelaxationMethod(F,dF,x1,size);
    
    InverseMatrix(dF, u, size);

    MatrixMultVector(u, F, buf, size);

    x2[0] = x1[0] - buf[0];
    x2[1] = x1[1] - buf[1];

    norm = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]));
    std::swap(x1, x2);
    iters++;

    while (norm > eps) {
        //std::cout << "Entered: " << std::endl;
        F[0] = f1(x1[0], x1[1]);
        F[1] = f2(x1[0], x1[1]);

        dF[0][0] = PartialDerivativeX(f1, x1[0], x1[1]);
        dF[0][1] = PartialDerivativeY(f1, x1[0], x1[1]);

        dF[1][0] = PartialDerivativeX(f2, x1[0], x1[1]);
        dF[1][1] = PartialDerivativeY(f2, x1[0], x1[1]);

        InverseMatrix(dF, u, size);

        MatrixMultVector(u, F, buf, size);

        x2[0] = x1[0] - buf[0];
        x2[1] = x1[1] - buf[1];

        norm = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]));
        std::swap(x1, x2);
        iters++;
        if (iters > 30) {
            std::cout << "La catastrophe!" << std::endl;
            break;
        }
    }
    return iters;
}


void NewtonNonlinearSearch(double*& mesh3DX, double*& mesh3DY, int& meshSize3DX, int& meshSize3DY,double*& F, double**& dF, double**& u, double*& buf, double*& x0, double*& x1, double*& x2,
    double (*f1)(double, double), double (*f2)(double, double), int size) {
    int t;
    std::vector<double> solsNX;
    std::vector<double> solsNY;
    std::vector<int> iters;

    std::ofstream fileOutput_1_10;
    fileOutput_1_10.open("res/1_10.txt");
    std::ofstream fileOutput_11_20;
    fileOutput_11_20.open("res/11_20.txt");
    std::ofstream fileOutput_21_30;
    fileOutput_21_30.open("res/21_30.txt");
    std::ofstream fileOutput_31_40;
    fileOutput_31_40.open("res/31_40.txt");
    std::ofstream fileOutputNoSol;
    fileOutputNoSol.open("res/noSol.txt");

    for (int i = 0; i < meshSize3DX; i++) {
        for (int j = 0; j < meshSize3DY - 1; j++) {
            t = HasSolution3D(mesh3DX[i], mesh3DY[j], mesh3DX[i], mesh3DY[j + 1], f1);
            switch (t) {
            case 0:
                solsNX.push_back(mesh3DX[i]);
                solsNY.push_back(mesh3DY[j]);
                //std::cout << "Solution 0 in "<< "["<< mesh3DX[j] <<","<< mesh3DY[j]<<"] " << "[" << mesh3DX[j] <<", "<< mesh3DY[j + 1]<<"]"  << std::endl;
                break;
            case 1:
                x0[0] = mesh3DX[i];
                x0[1] = (mesh3DY[j] + mesh3DY[j + 1]) / 2;
                iters.push_back(NewtonSystem(F, dF, u, buf, x0, x1, x2, f1, f2, size));
                solsNX.push_back(x1[0]);
                solsNY.push_back(x1[1]);
                DistributePointsToFiles(fileOutput_1_10, fileOutput_11_20, fileOutput_21_30, fileOutput_31_40, fileOutputNoSol, mesh3DX[i], (mesh3DY[j] + mesh3DY[j + 1]) / 2,iters[iters.size()-1]);
                //std::cout << "x0 = " << x0[0] << "y0 = " << x0[1] << " Solution " << "x = "<< x1[0]<< " y = "<< x1[1] <<", in " << "[" << mesh3DX[j] << "," << mesh3DY[j] << "] " << "[" << mesh3DX[j] << ", " << mesh3DY[j + 1] << "]" << std::endl;
                break;
            case 2:
                //fileOutputNoSol << mesh3DX[i] << " " << (mesh3DY[j] + mesh3DY[j + 1]) / 2 << "\n";
                //std::cout << "No solution in " << "[" << mesh3DX[j] << "," << mesh3DY[j] << "] " << "[" << mesh3DX[j] << ", " << mesh3DY[j + 1] << "]" << std::endl;
                break;
            }
        }
    }
    
    CheckForEqualElements3D(solsNX, solsNY,iters);

    std::cout << "Solution system(non unique)" << std::endl;
    for (int i = 0; i < solsNX.size(); i++) {
        std::cout << solsNX[i] << " " << solsNY[i] << std::endl;
    }

    fileOutput_1_10.close();
    fileOutput_11_20.close();
    fileOutput_21_30.close();
    fileOutput_31_40.close();
    fileOutputNoSol.close();
}

void RelaxationMethod(double*& F, double**& dF, double*& x, int size) {
    double tau = 0.01;
    double eps = 1e-6;
    double E[2][2]{};
    double x1[2]{};
    double x2[2]{};
    E[0][0] = 1.0;
    E[1][1] = 1.0;
    double sNorm;

    double S[2][2]{};
    bool flag = true;

    while (flag) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                S[i][j] = E[i][j] - dF[i][j] * tau;
            }
        }
        sNorm = CubicMatrixNorm(reinterpret_cast<double*>(S), size);
        std::cout << "sNorm " <<sNorm<< std::endl;
        if ( sNorm < 1) {
            break;
        }else {
            tau -= 0.001;
        }
        if (tau < 0.01) {
            //std::cout << "No right tau" << std::endl;
            flag = false;
            break;
        }
    }

    if (flag == false) {
        return;
    }

    x1[0] = x[0];
    x2[0] = x[0];

    x2[0] = x1[0] - tau * F[0];
    x2[1] = x1[1] - tau * F[1];

    double norm = sqrt((x2[0]-x1[0])*(x2[0] - x1[0]) + (x2[1] - x1[1])*(x2[1] - x1[1]));
    
    while (norm > eps) {
        x2[0] = x1[0] - tau * F[0];
        x2[1] = x1[1] - tau * F[1];
        norm = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]));
        std::swap(x1,x2);
    }

    x[0] = x1[0];
    x[1] = x1[1];
}