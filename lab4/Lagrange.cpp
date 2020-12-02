#include "Header.hpp"

double Interpolate(double x, Grid &grid)
{
    double ratio = 1;
    double res = 0;

    for(int i = 0; i < grid.length; ++i)
    {
        ratio = 1;
        for(int j = 0; j < grid.length; ++j)
        {
            if(i != j)
            {
                ratio *= (x - grid.points[j].x) / (grid.points[i].x - grid.points[j].x);
            }
        }
        res += grid.points[i].y * ratio;
    }

    return res;
}

void LagrangeInterpolate(Grid &initGrid, Grid &lagrange)
{
    for(int i = 0; i < lagrange.length; ++i)
    {
        lagrange.points[i].y = Interpolate(lagrange.points[i].x, initGrid);
    }
}

void LagrangeOut(Grid &lagrangeGrid, std::string fileName)
{
    std::ofstream fileOutput;
	fileOutput.open(fileName);

    for(int i = 0; i < lagrangeGrid.length; ++i)
    {
        fileOutput << lagrangeGrid.points[i].x << " " << lagrangeGrid.points[i].y << "\n";
    }

    fileOutput.close();
}

double CountError(Grid &grid, double (*f)(double))
{
    double err = 0;
    double iterErr = 0;
    for(int i = 0; i < grid.length; ++i)
    {
        iterErr = fabs(grid.points[i].y - f(grid.points[i].x));
        if(iterErr > err)
        {
            err = iterErr;
        }
    }

    return err;
}