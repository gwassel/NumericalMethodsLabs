#include "Header.hpp"


void LagrangeInterpolate(Grid &grid, Basis &basis, Polynomial &pLagrange, Polynomial &pBuffer)
{
    double a = 0, b = 0;// coeffitients of aX + b polynomial
    double denominator = 0;

    for(int k = 0; k < basis.length; ++k)
    {
        for(int j = 0; j < pLagrange.length; ++j)
        {
            if(k != j)
            {
                denominator = 1 / (grid.points[k].x - grid.points[j].x);
                a = denominator;
                b = grid.points[j].x * denominator;

                //multiplication of polynomials a*x + b
                for(int i = 0; i < j; ++i)
                {
                    pBuffer.coefficents[i] = basis.polynomials[k].coefficents[i] * b;
                    pBuffer.coefficents[i + 1] = basis.polynomials[k].coefficents[i + 1] * a;
                }

                for(int i = 0; i < j + 1; ++i)
                {
                    basis.polynomials[k].coefficents[i] = pBuffer.coefficents[i];
                }
            }
        }
        for(int j = 0; j < pLagrange.length; ++j)
        {
            basis.polynomials[k].coefficents[j] *= grid.points[k].y;
        }
    }
    
    for(int i = 0; i < basis.length; ++i)
    {
    //sum
    }


}

void SplineInterpolate()
{
    
}
 
