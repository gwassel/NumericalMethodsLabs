#include "Header.h"

double Derivative(double (*f)(double), double x) {
    double epsilon = 1e-6;
    return (f(x + epsilon) - f(x)) / epsilon;
}

double SecondDerivative(double (*f)(double), double x) {
    double epsilon = 1e-6;
    return (f(x + epsilon) - 2 * f(x) + f(x - epsilon)) / (epsilon * epsilon);
}

double PartialDerivativeX(double (*f)(double, double), double x, double y) {
    double epsilon = 1e-6;
    return (f(x + epsilon, y) - f(x, y)) / epsilon;
}

double PartialDerivativeY(double (*f)(double, double), double x, double y) {
    double epsilon = 1e-6;
    return (f(x, y + epsilon) - f(x, y)) / epsilon;
}

// 0 begin of otrezok koren
// 1 koren na otrezke
// 2 kornya net
int HasSolution2D(double xi, double xj, double (*f)(double)) {
    double epsilon = 1e-6;
    if (fabs(f(xi)) < epsilon) return 0;
    if (f(xi) * f(xj) < 0) {
        return 1;
    }
    return 2;
}

// 0 begin of otrezok koren
// 1 koren na otrezke
// 2 kornya net
//возможно требуется модификация
int HasSolution3D(double xi, double yi, double xj, double yj, double (*f)(double,double)) {
    double epsilon = 1e-6;
    if (fabs(f(xi,yi)) < epsilon) return 0;
    if (f(xi,yi) * f(xj,yj) < 0) {
        return 1;
    }
    return 2;
}


void MakeMesh2D(double A, double B, double*& table, int& size) {
    double step = (B - A) / (size - 1.0);
    for (int i = 0; i < size; i++) {
        table[i] = A + i * step;
    }
}

void PrintTable(double*& table, int& size) {
    std::cout << "Mesh: " << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << table[i] << " ";
    }
    std::cout << std::endl;
}

void CheckForEqualElements(std::vector<double>& sols) {
    double eps = 1e-6;
    bool flag = true;
    while (flag) {
        flag = false;
        for (int i = 0; i < sols.size(); i++) {
            for (int j = 0; j < sols.size(); j++) {
                if (i == j) continue;
                if (fabs(sols[i] - sols[j]) < eps) {
                    sols.erase(sols.begin() + j);
                    flag = true;
                    break;
                }
            }
        }
    }
}

void CheckForEqualElements3D(std::vector<double>& solsNX, std::vector<double>& solsNY, std::vector<int>& iters) {
    double eps1 = 1e-8;
    bool flag = true;
    while (flag) {
        flag = false;
        for (int i = 0; i < solsNX.size(); i++) {
            for (int j = 0; j < solsNX.size(); j++) {
                if (i == j) continue;
                if (fabs(solsNX[i] - solsNX[j]) < eps1 && fabs(solsNY[i] - solsNY[j]) < eps1) {
                    solsNX.erase(solsNX.begin() + j);
                    solsNY.erase(solsNY.begin() + j);
                    iters.erase(iters.begin() + j);
                    flag = true;
                    break;
                }
            }
        }
    }
}


void DistributePointsToFiles(std::ofstream& f1_10, std::ofstream& f11_20, std::ofstream& f21_30, std::ofstream& f31_40, std::ofstream& fnoSols,
    double x, double y, int iters) {
    //std::cout << "Iters: " << iters << std::endl;
    /*if (iters >= 1 && iters < 11) {
       f1_10 << x << " " << y << "\n";
    }
    else if(iters >= 11 && iters < 21){
            f11_20 << x << " " << y << "\n";
    }
    else if (iters>= 21 && iters < 31) {
            f21_30 << x << " " << y << "\n";
    }
    else if (iters>= 31 && iters < 41) {
            f31_40 << x << " " << y << "\n";*/
    if (iters >= 1 && iters < 3) {
        f1_10 << x << " " << y << "\n";
    }
    else if (iters >= 3 && iters < 5) {
        f11_20 << x << " " << y << "\n";
    }
    else if (iters >= 5 && iters < 7) {
        f21_30 << x << " " << y << "\n";
    }
    else if (iters >= 7 && iters < 10) {
        f31_40 << x << " " << y << "\n";
    }
    else if (iters >= 10) {
        fnoSols << x << " " << y << "\n";
    }
}

