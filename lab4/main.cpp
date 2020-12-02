#include "Header.hpp"

int main()
{
    //init variables
    size_t n = 7;
    size_t size = n + 1;
    size_t testSize = 1024 + 1;

    double leftBorder = -1;
    double rightBorder = 1;

    //non init variables
    Grid uniformGrid = Grid(size);
    Grid ChebishevGrid = Grid(size);

    Grid lagrangeGrid1 = Grid(testSize);
    Grid lagrangeGrid2 = Grid(testSize);

    Grid testGrid = Grid(testSize);

    // Grid splineGrid = Grid(6);
    // splineGrid.points[0] = {1,1.0002};
    // splineGrid.points[1] = {2,1.0341};
    // splineGrid.points[2] = {3,0.6};
    // splineGrid.points[3] = {4,0.40105};
    // splineGrid.points[4] = {5,0.1};
    // splineGrid.points[5] = {6,0.23975};

    // Grid splineGrid = Grid(6);
    // splineGrid.points[0] = {0,1};
    // splineGrid.points[1] = {3,2};
    // splineGrid.points[2] = {8,3};
    // splineGrid.points[3] = {15,4};
    // splineGrid.points[4] = {24,5};
    // splineGrid.points[5] = {35,6};

    // Grid splineGrid = Grid(3);
    // splineGrid.points[0] = {0, 1};
    // splineGrid.points[1] = {0.6, cos(0.6 * 0.6)};
    // splineGrid.points[2] = {0.9, cos(0.9 * 0.9)};
    // splineGrid.points[3] = {1.2, cos(1.2 * 1.2)};


    Spline s1;


    double (*f)(double) = f7; //init function to interpolate

    MakeMesh(leftBorder, rightBorder, lagrangeGrid1, 1, f0);
    MakeMesh(leftBorder, rightBorder, lagrangeGrid2, 1, f0);

    MakeMesh(leftBorder, rightBorder, uniformGrid, 1, f);
    MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2, f);

    LagrangeInterpolate(uniformGrid, lagrangeGrid1);
    LagrangeInterpolate(ChebishevGrid, lagrangeGrid2);

    LagrangeOut(lagrangeGrid1, "lagrange1");
    LagrangeOut(lagrangeGrid2, "lagrange2");


    double err1 = 0;
    double err2 = 0;
    err1 = CountError(lagrangeGrid1, f);
    err2 = CountError(lagrangeGrid2, f);

    std::cout << "err on uniform: " << err1 << "\nerr on Chebishev: " << err2 << "\n";
//lagrange ends


    s1 = Spline(uniformGrid);
    TridiagonalMatrix splineMatrix(s1.h, s1.g, s1.n + 1);
    splineMatrix.eval();
    splineMatrix.run();
    s1.RecountCoefficents(splineMatrix);
    // for(int i = 0; i < s1.n; ++i)
    //     std::cout << "a: " << s1.a[i] << " b: " << s1.b[i] << " c: " << s1.c[i] << " d: " << s1.d[i] << "\n";

    WriteSpline("splineout", s1);
    double errspline = CountError(s1, f);
    std::cout << "err on spline: " << errspline << "\n";
    // std::cout << s1.eval(0) << "\n";
    return 0;
}
