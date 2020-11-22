#include "Header.hpp"


bool comparator(Point &p1, Point &p2)
{
	return p1.x < p2.y;
}

//uniform mesh - 1, Chebyshev mesh - 2
//PointVector  size = n + 1, must count in x0!
void MakeMesh(double x0, double xN, Grid &grid, int MeshType, double (*f)(double)) {
	double step;
	const double PI = 3.1415;
	if (grid.length != 0) {
		step = (xN - x0) / (grid.length - 1.0);
	}
	else {
		return;
	}

	switch (MeshType) {
	case 1:
		grid.points[0] = x0;
		for (int i = 1; i < grid.length; i++) {
			grid.points[i].x = x0 + i * step;
		}
		break;
	case 2:
		double A = (xN + x0) / 2.0;
		double B = (xN - x0) / 2.0;
		double C = 2.0 * (grid.length + 1.0);
		std::cout << grid.length << std::endl;
		for (int i = 0; i < grid.length; i++) {
			grid.points[i].x = A + B * cos(PI * (2.0 * i + 1.0) / C);
		}
		//sort
		std::sort(grid.points, grid.points + grid.length, comparator);
	}

	for (int i = 0; i < grid.length; i++) {
		grid.points[i].y = f(grid.points[i].x);
	}

	WriteCoords("res/coord.txt", grid);
}
