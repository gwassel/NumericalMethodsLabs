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
            }
        }
    }
}

void SplineInterpolate()
{
    
}
 
