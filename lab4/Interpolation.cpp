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

void LagrangeInterpolate(Grid &grid, Basis &basis, Polynomial &pLagrange, Polynomial &pBuffer, Polynomial &monomial, std::string label /*= ""*/)
{
    std::cout << label << "\n";
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
        pLagrange.coefficents[i] = 0;

        for(int j = 0; j < basis.length; ++j)
        {
            pLagrange.coefficents[i] += basis.polynomials[j].coefficents[i] * grid.points[j].y;
        }
    }
    for(int i = 0; i < pLagrange.length; ++i)
        std::cout << pLagrange.coefficents[i] << "x^" << i << ", ";
    std::cout << "\n";
}

void SplineInterpolate()
{
    
}
 
double CheckOnGridPoints(Polynomial &p1, Grid &grid)
{
    double err = 0;
    for(int i = 0; i < p1.length; ++i)
    {
        err += fabs(p1.eval(grid.points[i].x) - grid.points[i].y);
    }
    return err;
}