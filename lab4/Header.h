#pragma once
#include <string>
void MakeMesh(double x0, double xN, double*& PointVector, int size, int MeshType);
void MergeSort(double*& A, int first, int last, int size);
void Merge(double*& A, int first, int last, int size);
void WriteCoords(const std::string fileNameOutput, double*& vectorX, double*& vectorY, const size_t n);