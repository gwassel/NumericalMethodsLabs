#include "Header.h"
double f1(double x)
{
    return std::tan(x);
}

double f2(double x)
{
    double r1 = tan(5 * x * x + 3 * x - 2);
    double r2 = exp((x * x * x + 6 * x * x + 12 * x + 8) / (2 * x * x + 8 * x + 7));
    return r1 + r2 - 2.0;
}

double f3(double x, double y) {
    return x * x - 1;
}

double f4(double x, double y) {
    return x * x * x - x * x - 4 + y;
}

double ftest(double x, double y) {
    return sin(x * y) * cos(x);
}

//System 1
double f5(double x, double y) {
    return x * x - y * y - 15;
}
double f6(double x, double y) {
    return x * y + 4;
}

//System 2
double f7(double x, double y) {
    return x * x + y * y + x + y - 8;
}
double f8(double x, double y) {
    return x * x + y * y + x * y - 7;
}

//System var4
double f9(double x, double y) {
    return 3 * x * x - 2 * y + 8;
}
double f10(double x, double y) {
    return x + y - 4;
}

//System var1
double f11(double x, double y) {
    return 2 * (x + y) * (x + y) + (x - y) * (x - y) - 8;
}
double f12(double x, double y) {
    return 5 * x * x + (y - 3) * (y - 3) - 9;
}

//One eq, ¹2
double f13(double x)
{
    return sqrt(x+1) - 1;
}

//One eq, ¹3
double f14(double x)
{
    return 35 * x * x * x - 67 * x * x - 3 * x + 3;
}

double f15(double x)
{
    return x;
}
