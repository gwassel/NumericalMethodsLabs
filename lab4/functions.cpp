#include "Header.hpp"

double f0(double x)
{
    return -1;
}

double f1(double x)
{
    return sin(x);
}

double f2(double x)
{
    return 1.0 / (1 + x * x);
}

double f3(double x)
{
    return 1;
}

double f4(double x)
{
    return x * x;
}

double f5(double x)
{
    double r1 = tan( (pow(x, 4) + 2 * pow(x, 2) - 2 * x + sqrt(2) + 1) / 8);
    double r2 = sin( (4 * pow(x, 3) - 7 * x - 9) / (20 * x + 28) );
    return r1 + r2;
}

double f6(double x)
{
    return 1.0 / (1 + 25 * x * x);
}

double f7(double x)
{
    return sin(3.1415926535 * x);
}