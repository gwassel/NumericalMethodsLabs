#include "Header.hpp"

int main()
{
    //init variables
    size_t n = 30;
    size_t size = n + 1;
    size_t testSize = 1024 + 1;

    double leftBorder = -1;
    double rightBorder = 1;

    //non init variables
    Grid uniformGrid = Grid(size);
    Grid ChebishevGrid = Grid(size);
    Grid testGrid = Grid(testSize);

    // Grid splineGrid = Grid(6);
    // splineGrid.points[0] = {1,1.0002};
    // splineGrid.points[1] = {2,1.0341};
    // splineGrid.points[2] = {3,0.6};
    // splineGrid.points[3] = {4,0.40105};
    // splineGrid.points[4] = {5,0.1};
    // splineGrid.points[5] = {6,0.23975};

    Grid splineGrid = Grid(3);
    splineGrid.points[0] = {0,1};
    splineGrid.points[1] = {3,2};
    splineGrid.points[2] = {8,3};

    Polynomial p1 = Polynomial(size);
    Polynomial p2 = Polynomial(size);
    Polynomial pBuffer = Polynomial(size);
    Polynomial monomial = Polynomial(2);
    
    Basis basis1 = Basis(size);
    Basis basis2 = Basis(size);

    Spline s1;


    double (*f)(double) = f1; //init function to interpolate

    // MakeMesh(leftBorder, rightBorder, uniformGrid, 1, f);
    // MakeMesh(leftBorder, rightBorder, uniformGrid, 1, f);
    MakeMesh(leftBorder, rightBorder, ChebishevGrid, 2, f);

    ReadInit();
    ReadData();

    // LagrangeInterpolate(uniformGrid, basis1, p1, pBuffer, monomial, "Lagrange Polynomial on uniform grid");
    LagrangeInterpolate(ChebishevGrid, basis2, p2, pBuffer, monomial, "Lagrange Polynomial on Chebishev grid");

    s1 = Spline(splineGrid);
    TridiagonalMatrix splineMatrix(s1.h, s1.g, s1.n + 1);
    splineMatrix.eval();
    splineMatrix.run();
    s1.RecountCoefficents(splineMatrix);
    for(int i = 0; i < s1.n; ++i)
        std::cout << "a: " << s1.a[i] << " b: " << s1.b[i] << " c: " << s1.c[i] << " d: " << s1.d[i] << "\n";


    SplineInterpolate();

    WriteData();

    // test(p1, uniformGrid, "Lagrange on uniform grid");
    test(p2, ChebishevGrid, "Lagrange on Chebishev grid");


    // double err;
    // err = CountError(p1, testGrid, f, leftBorder, rightBorder);
    // std::cout << "err on " << testSize - 1 << " grid = " << err << "\n";
    // err = CountError(p2, testGrid, f, leftBorder, rightBorder);
    // std::cout << "err on " << testSize - 1 << " grid = " << err << "\n";
    
    // WriteSpline("spline.xdd", s1);
    // WritePolynomial("lagrangexd", p2, testGrid);



    // std::cout << s1.eval(0) << "\n";
    return 0;
}
