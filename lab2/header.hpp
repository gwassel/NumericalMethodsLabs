#include <iostream>
#include <fstream>


//matrix functions
void MatrixMult(double**&, double**&, double**&, const size_t);
void MatrixMult(double**&, double*&, double*&, const size_t);

void VectorCopy(double*&, double*&, const size_t);
void VectorAdd(double*&, double*&, double*&, const size_t);

void MatrixEDiff(double**&, double**&, const size_t);


//Molochkov functions
bool IfStop(double**&, double*&, double*&, double*&, const double, const size_t);
void Iterations(double**&, double*&, double*&, double*&, double*&, const double, const size_t);
void GetCY(double**&, double*&, double, double**&, double*&, const size_t);
