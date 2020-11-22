#include "Header.hpp"


void PolynomialMultiplication(Polynomial &p1, Polynomial &p2, Polynomial &pResult)
{
    for(int i = 0; i < pResult.length; ++i)
    {
        pResult.coefficents[i] = 0;
    }

    for(int i = 0; i < p1.length; ++i)
    {
        for(int j = 0; j < p2.length; ++j)
        {
            pResult.coefficents[i + j] += p1.coefficents[i] * p2.coefficents[j];
        }
    }
}

void LagrangeInterpolate(Grid &grid, Basis &basis, Polynomial &pLagrange, Polynomial &pBuffer, Polynomial &monomial)
{
    double a = 0, b = 0; // coeffitients of aX + b polynomial
    double denominator = 0;

    for(int k = 0; k < basis.length; ++k)
    {
        for(int j = 0; j < pLagrange.length; ++j)
        {
            if(j != k)
            {
                denominator = 1 / (grid.points[k].x - grid.points[j].x);

                monomial.coefficents[1] = denominator;
                monomial.coefficents[0] = - grid.points[j].x * denominator;

                PolynomialMultiplication(basis.polynomials[k], monomial, pBuffer);

                for(int i = 0; i < pBuffer.length; ++i)
                {
                    basis.polynomials[k].coefficents[i] = pBuffer.coefficents[i];
                }
            }
        }
    }
    for(int i = 0; i < pLagrange.length; ++i)
    {
        a = 0;

        for(int j = 0; j < basis.length; ++j)
        {
            a += basis.polynomials[j].coefficents[i] * grid.points[j].y;
        }

        pLagrange.coefficents[i] = a;
    }
}

void SplineInterpolate()
{
    
}
 
double test(Polynomial &p1, Grid &grid)
{
    double err = 0;
    for(int i = 0; i < p1.length; ++i)
    {
        err += fabs(p1.eval(grid.points[i].x) - grid.points[i].y);
    }
    return err;
}