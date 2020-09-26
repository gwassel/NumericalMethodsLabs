#include "header.hpp"

void WriteMatrixTeX(std::string fileOutputName, const double** &matrix, const size_t n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileOutputName);

    fileOutput << "begin{matrix}\n";

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            fileOutput << matrix[i][j];
            if (j != n - 1)
                fileOutput << "& ";
        }
        fileOutput << "\\ \\ \n";
    }
}
