#include "Header.h"
//сделать проверку на точки разрыва
//sin^2(x) не находит sin(x*x)
//дописатьм етод релакцации
int main()
{
    //Section
    double leftBorder = -1.0;
    double rigthBorder = 2.0;

    double leftBorderX = -10;
    double rigthBorderX = 10;
    double leftBorderY = -10;
    double rigthBorderY = 10;

    //Function for bisection and Newton method
    double (*f)(double) = f15;

    //Size of table(number of points) == number of sections + 1
    int meshSize2D = 17;
    int meshSize3DX = 100;
    int meshSize3DY = 100;

    //Mesh
    double* mesh2D;
    double* mesh3DX;
    double* mesh3DY;

    //Vector function
    double* F;

    //Size of vector function
    int sizeF = 2;

    //Inversed F
    double** inversedF;
    double** dF;

    //Vectors for Newton system method
    double* buf;
    double* x0;
    double* x1;
    double* x2;

    AllocateMemory(mesh2D, mesh3DX, mesh3DY, F, inversedF, dF, buf, x0, x1, x2,
        sizeF, meshSize2D, meshSize3DX, meshSize3DY);

    MakeMesh2D(leftBorder, rigthBorder, mesh2D, meshSize2D);

    MakeMesh2D(leftBorderX, rigthBorderX, mesh3DX, meshSize3DX);
    MakeMesh2D(leftBorderY, rigthBorderY, mesh3DY, meshSize3DY);

    PrintTable(mesh2D, meshSize2D);

    CalculateOneEquation(mesh2D, meshSize2D, leftBorder, rigthBorder, f);
    NewtonNonlinearSearch(mesh3DX, mesh3DY, meshSize3DX, meshSize3DY, F, dF, inversedF, buf, x0, x1, x2, f9, f10, sizeF);
    //RelaxationMethod(F,dF,sizeF);
}
