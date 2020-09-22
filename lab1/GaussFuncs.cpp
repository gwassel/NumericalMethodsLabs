#include "funcs.hpp"

const double epsilon = 1e-12;

int EMain(const int Esize, my_type** EA, my_type** EIA, my_type* EB, my_type* EBF, my_type* EdB, my_type* EX1, 
        my_type* EX2, my_type Econd, my_type Edx, my_type Edb, bool EnoProblem, std::tuple<my_type, my_type, my_type> Ed11,
        std::tuple<my_type, my_type> Ec11, const std::string EpathConfig, const std::string EpathData, std::string EfileNameMatrix,
        std::string EfileNameVector, const std::string EfileMatrixAName, const std::string EfileMatrixAIName, const std::string EfileVectorX1Name,
        const std::string EfileVectorX2Name, const std::string EfileVectorBName, const std::string EfileParamsName, std::string EfolderName)
{
    //1--------------------------------------------------------------------------------
    EAllocateMemory(EA,EIA,EB,EdB,EX1, EX2,EBF,Esize);
    EfileNameMatrix = EpathData + "matrix1";
    EfileNameVector = EpathData + "vector11";
    EReadData(EfileNameMatrix, EfileNameVector,EA,EB,Esize);
 
    EfolderName = "results/res" + std::to_string(1);
    ECalculations(EfolderName + EfileMatrixAName, EfolderName + EfileMatrixAIName, EfolderName + EfileVectorBName,
        EfolderName + EfileVectorX1Name, EfolderName + EfileVectorX2Name, EfolderName + EfileParamsName, EA, EIA, EB, EdB, EX1,EX2,EBF, Esize);
 
    //3--------------------------------------------------------------------------------
    EfileNameMatrix = EpathData + "matrix2";
    EfileNameVector = EpathData + "vector21";
    EfolderName = "results/res" + std::to_string(3);
    EReadData(EfileNameMatrix, EfileNameVector, EA, EB, Esize);
 
    ECalculations(EfolderName + EfileMatrixAName, EfolderName + EfileMatrixAIName, EfolderName + EfileVectorBName,
        EfolderName + EfileVectorX1Name, EfolderName + EfileVectorX2Name, EfolderName + EfileParamsName, EA, EIA, EB, EdB, EX1,EX2,EBF, Esize);
 
    EFreeMemory(EA, EIA, EB, EdB, EX1, EX2, EBF, Esize);

    return 0;
}

int EAllocateMemory(my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, my_type*& BF, const int& n)
{
    A = new my_type * [n];
    BF = new my_type[n];

    for (int i = 0; i < n; ++i){
        A[i] = new my_type[n];
    }

    IA = new my_type * [n];

    for (int i = 0; i < n; ++i)
    {
        IA[i] = new my_type[n]{};
    }

    B = new my_type[n];
    X1 = new my_type[n];
    X2 = new my_type[n];

    dB = new my_type[n];
    return 0;
}


int EReadData(const std::string fileNameMatrix, const std::string fileNameVector, my_type**& matrixA, my_type*& vectorB, const int& n)
{
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with matrix is not open\n";
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixFile >> matrixA[i][j];
        }
    }

    matrixFile.close();

    std::ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if (!vectorFile.is_open())
    {
        std::cerr << "Error: file with vector is not open\n";
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        vectorFile >> vectorB[i];
    }

    return 0;
}

void EprintEquation(my_type**& A, my_type*& B, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << B[i];
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void EprintMatrix(my_type**& A, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void EprintVector(my_type*& B, const int& size, std::string s){
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << B[i] << std::endl;
    }
    std::cout << std::endl;
}

void EGaussMethod(my_type**& A1, my_type*& B, my_type*& X, const int& size, bool& flag) {
    my_type m;
    my_type k;
    bool noProblems = true;
    my_type dub;
    std::vector<std::tuple<int, int>> permutations = {};
    my_type** A = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        A[i] = new my_type[size];
    }

    for (int i = 0; i < size; ++i) {
        X[i] = B[i];
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++ ) {
            A[i][j] = A1[i][j];
        }
    }

    /*for (int i = 0; i < size - 1; i++) {
        if (MatrixIsPrepared(A,B,i,permutations,size)) {
            m = A[i][i];
            B[i] /= m;
            for (int j = i; j < size; j++) {
                for (int i1 = i + 1; i1 < size; i1++) {
                    k = A[i1][i];
                    B[i1] -= k * B[i];
                    for (int j1 = i; j1 < size; j1++) {
                        A[i1][j1] -= k * (A[i][j1] / m);
                    }
                }
                A[i][j] /= m;
            }
           }*/
    for (int i = 0; i < size - 1; i++) {
        if (EMatrixIsPrepared(A, B, i, permutations, size)) {
            m = A[i][i];
            B[i] /= m;

            for (int j = i; j < size; j++) {
                A[i][j] /= m;
            }

            for (int i1 = i + 1; i1 < size; i1++) {

                k = A[i1][i];
                B[i1] -= k * B[i];

                for (int j1 = i; j1 < size; j1++) {
                    A[i1][j1] -= k * A[i][j1];
                }
            }
        }
        else {
            noProblems = false;
            break;
        }
    }
    if (EMatrixIsPrepared(A,B,size - 1, permutations, size) == false) {
        noProblems = false;
    }
    if (noProblems == true) {
        B[size - 1] /= A[size - 1][size - 1];
        A[size - 1][size - 1] /= A[size - 1][size - 1];
        for (int i = size - 1; i > -1; i--) {
            for (int j = i-1; j > -1; j--) {
                        B[j] -= A[j][i] * B[i];
                        A[j][i] = 0.0;//A[j][i] -= A[j][i];
            }
        }
        /*B[size - 1] /= A[size - 1][size - 1];
        A[size - 1][size - 1] /= A[size - 1][size - 1];
        for (int i = size - 1; i > -1; i--) {
           m = A[i][i];
            B[i] /= A[i][i];
        for (int j = i - 1; j > -1; j--) {
            /*for (int i1 = i - 1; i1 > -1; i1--) {
                k = A[i1][i];
                B[i1] -= k * B[i];
                for (int j1 = i; j1 > -1; j1--) {
                //A[i1][j1] -= A[i][j1] / m;
            //
     //   }
    //}
    //A[i][j] /= m;
        }
    }*/
        EDiagonalizeEquation(A,B,X,size,permutations);
        for (int i = 0; i < size; i++) {
            dub = X[i];
            X[i] = B[i];
            B[i] = dub;
        }
    }
    flag = noProblems;
    permutations.clear();
    for (int i = 0; i < size; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
}

bool EMatrixIsPrepared(my_type**& A, my_type*& B, const int& i, std::vector<std::tuple<int, int>>& permutations, const int& size) {
    std::tuple<int, int> t;
    t = ESearchMax(A,i,size);

    if (std::get<0>(t) != i) {
        ESwapRows(A,B,i, std::get<0>(t));
    }
    //if (std::get<1>(t) != i) {
    //    ESwapColomns(A,i, std::get<1>(t),size);
    //    permutations.push_back(std::make_tuple(i, std::get<1>(t)));
    //}
    if (EisDegenerate(A,i,size) == false) {
        return true;
    }
    else return false;
}

void EDiagonalizeEquation(my_type**& A, my_type*& B, my_type*& X, const int& size, std::vector<std::tuple<int, int>>& permutations) {
    for (int i = permutations.size() - 1; i > -1; i--) {
        ESwapColomns(A, std::get<0>(permutations[i]), std::get<1>(permutations[i]), size);
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (A[i][j] > epsilon && i != j) { 
                ESwapRows(A, B, i, j);
            }
        }
    }
    //permutations.clear();
}

std::tuple<int, int> ESearchMax(my_type**& A, const int& currentMinor, const int& size) {
    /*int max_i = currentMinor;
    int max_j = currentMinor;
    my_type max = fabs(A[max_i][max_j]);
    for (int i = currentMinor; i < size; i++) {
        for (int j = currentMinor; j < size; j++) {
            if (fabs(A[i][j]) > max) {
                max_i = i;
                max_j = j;
                max = fabs(A[max_i][max_j]);
            }
        }
    }
    return std::make_tuple(max_i, max_j);*/

    int max_i = currentMinor;
    my_type max = fabs(A[max_i][currentMinor]);
    for (int i = currentMinor; i < size; i++) {
            if (fabs(A[i][currentMinor]) > max) {
                max_i = i;
                max = fabs(A[max_i][currentMinor]);
        }
    }
    return std::make_tuple(max_i, currentMinor);
}

void ESwapRows(my_type**& A, my_type*& B, const int& i1, const int& i2) {
    //my_type* dub1 = new my_type;
    std::swap(A[i1], A[i2]);
    //dub1 = A[i1];
    //A[i1] = A[i2];
    //A[i2] = dub1;

    std::swap(B[i1], B[i2]);
    //my_type dub2 = B[i1];
    //B[i1] = B[i2];
    //B[i2] = dub2;
    //std::cout << "SwapRows " << i1 << " " << i2 << std::endl;
    //printEquation("Swap Rows!");
}

void ESwapColomns(my_type**& A, const int& j1, const int& j2, const int& size) {
    my_type dub;
    for (int i = 0; i < size; i++) {        
        dub = A[i][j1];
        A[i][j1] = A[i][j2];
        A[i][j2] = dub;
    }

    //std::cout << "SwapColomns "<<j1<<" "<<j2 << std::endl;
    //printEquation("Swap Colomns!");
}

bool EisDegenerate(my_type**& A, const int& i,const int& size) {
    int k = 0;
    for (int j = 0; j < size; j++) {
        if (fabs(A[i][j]) < epsilon) {
            k++;
        }
    }

    /*if (fabs(A[i][i])<epsilon) {
        return true;
    }*/

    if (k == size) {
        return true;
    }
    return false;
}

my_type ECubicVectorNorm(my_type* &p, const int& size) {
    my_type sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(p[i]);
    }
    return sum;
}

my_type ECubicMatrixNorm(my_type** &p, const int& size) {
    my_type sum;
    my_type maxSum = 0.0;
    for (int j = 0; j < size; j++) {
        maxSum += fabs(p[0][j]);
    }
    for (int i = 1; i < size; i++) {
        sum = 0.0;
        for (int j = 0; j < size; j++) {
            sum += fabs(p[i][j]);
        }
        if (sum > maxSum)
            maxSum = sum;
    }
    return maxSum;
}

my_type EOctahedralVectorNorm(my_type* &p, const int& size) {
    my_type max = fabs(p[0]);
    for (int i = 1; i < size; i++) {
        if (max < fabs(p[i]))
            max = fabs(p[i]);
    }
    return max;
}

my_type EOctahedralMatrixNorm(my_type** &p, const int& size) {
    my_type sum;
    my_type maxSum = 0.0;
    for (int i = 0; i < size; i++) {
        maxSum += fabs(p[i][0]);
    }
    for (int i = 1; i < size; i++) {
        sum = 0.0;
        for (int j = 0; j < size; j++) {
            sum += fabs(p[j][i]);
        }
        if (sum > maxSum)
            maxSum = sum;
    }
    return maxSum;
}

std::tuple<my_type, my_type> EConditionNumber(my_type**& A, my_type**& IA, const int& size) {

    //my_type c;

    //c = CubicMatrixNorm(IA, size) * CubicMatrixNorm(A, size);
    //std::cout << "Condition number with cubic matrix norm:" << c << std::endl;

    //c = OctahedralMatrixNorm(IA, size) * OctahedralMatrixNorm(A, size);
    //std::cout << "Condition number with octahedral matrix norm:" << c << std::endl;
    //std::cout << std::endl;
    return std::make_tuple(ECubicMatrixNorm(IA, size) * ECubicMatrixNorm(A, size), EOctahedralMatrixNorm(IA, size) * EOctahedralMatrixNorm(A, size));
}

int EMatrixMult(my_type**& matrixA, my_type**& matrixB, const int& n)
{
    my_type** matrixResult = new my_type*[n];
    for (int i = 0; i < n; i++) {
        matrixResult[i] = new my_type[n]{};
    }
    my_type sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            sum = 0.0;
            for (int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (matrixResult[i][j] < epsilon) {
                matrixResult[i][j] = 0.0;
            }
        }
    }

    EprintMatrix(matrixResult, n, "A*A-1");

    for (int i = 0; i < n; ++i)
    {
        delete[] matrixResult[i];
    }
    delete[] matrixResult;

    return 0;
}//matrixA * matrixB


std::tuple<my_type,my_type,my_type> Ediscrepancy(my_type**& A, my_type*& B, my_type*& X, my_type*& BF, const int& size) {
    my_type* D = new my_type[size]{};
    my_type sum;

    for (int i = 0; i < size; i++)
    {
        sum = 0.0;
        for (int j = 0; j < size; j++)
        {
            sum += A[i][j] * X[j];
        }
        D[i] = sum;
    }
    //printVector(D,size,"D");
    //printVector(B,size,"B");
    //sum = 0.0;
    //for (int i = 0; i < size; i++) {
    //    sum += pow(D[i] - B[i], 2);
   // }

    //std::cout << "Discrepancy(euclid) " << sqrt(sum) << std::endl;
   // std::cout << "Discrepancy(cubic) " << CubicVectorNorm(D, size) << std::endl;
    //std::cout << "Discrepancy(octahedral) " << OctahedralVectorNorm(D, size) << std::endl;
   // std::cout << std::endl;
    std::tuple<my_type, my_type, my_type> t = std::make_tuple(sqrt(sum), ECubicVectorNorm(EminusVectors(D, B,BF, size), size), EOctahedralVectorNorm(EminusVectors(D, B,BF, size), size));
    delete[] D;
    return t;
}

int EFreeMemory(my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, my_type*& EBF, const int& n)
{
    delete[] B;
    delete[] dB;
    delete[] X1;
    delete[] X2;
    delete[] EBF;

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
        delete[] IA[i];
    }
    delete[] A; 
    delete[] IA;

    return 0;
}

my_type*& EminusVectors(my_type*& X1, my_type*& X2, my_type*& BF, const int& size) {
    for (int i = 0; i < size; i++) {
        BF[i] = X1[i] - X2[i];
    }
    return BF;
}

my_type EmarkCondition(my_type*& X1, my_type*& X2, my_type*& B1, my_type*& B2, my_type& dx, my_type& db, const int& size) {
    
    my_type normX1 = ECubicVectorNorm(X1,size);
    my_type normX2 = ECubicVectorNorm(X2, size);
    my_type normB1 = ECubicVectorNorm(B1,size);
    my_type normB2 = ECubicVectorNorm(B2,size);
    //my_type sum;
    /*
    sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(X1[i], 2);
    }
    normX1 = sqrt(sum);

    sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(X2[i], 2);
    }
    normX2 = sqrt(sum);

    sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(B1[i], 2);
    }
    normB1 = sqrt(sum);

    sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(B2[i], 2);
    }
    normB2 = sqrt(sum);
    */
    my_type dX = normX1 / normX2;
    dx = dX;
    my_type dB = normB1 / normB2;
    db = dB;
    return  dX / dB;
}

int EWriteData(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP,my_type**& matrixA, my_type**& matrixIA,
    my_type*& vectorX1, my_type*& vectorX2, my_type*& vectorB, const my_type dis1,
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx, const int& n)
{
    EWriteMatrix(fileNameA, "matrix A", matrixA, n);
    EWriteMatrix(fileNameIA, "matrix IA", matrixIA, n);

    EWriteVector(fileNameB, "vector B", vectorB, n);
    EWriteVector(fileNameX1, "vector X1", vectorX1, n);
    EWriteVector(fileNameX2, "vector X2", vectorX2, n);

    EWriteParameters(fileNameP, dis1, dis2, dis3, condb, condNum1, condNum2, db, dx);
    return 0;
}

int EWriteMatrix(const std::string fileNameOutput, const std::string label, my_type**& matrix, const int& n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    fileOutput << label << "\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j) {
            fileOutput << matrix[i][j] << " ";
        }
        fileOutput << "\n";
    }
    fileOutput << "\n";

    fileOutput.close();

    return 0;
}

int EWriteVector(std::string fileNameOutput, const std::string label, my_type*& vector, const int& n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    fileOutput << label << "\n";
    for (int i = 0; i < n; ++i)
    {
        fileOutput << vector[i] << "\n";
    }
    fileOutput << "\n";

    fileOutput.close();

    return 0;
}

int EWriteParameters(std::string fileNameOutput, const my_type dis1,
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx){
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    //fileOutput << "Discrepancy(euclid) " << dis1 << "\n";
    fileOutput << "Discrepancy(cubic) " << dis2 << "\n";
    fileOutput << "Discrepancy(octahedral) " << dis3 << "\n";
    fileOutput << "\n";

    fileOutput << "Condition number is bigger than " << condb << "\n";
    fileOutput << "\n";

    fileOutput << "Condition number with cubic matrix norm:" << condNum1 << "\n";
    fileOutput << "Condition number with octahedral matrix norm:" << condNum2 << "\n";
    fileOutput << "\n";

    fileOutput << "Relative error of X is  " << dx << "\n";
    fileOutput << "Relative error of B is " << db << "\n";
    fileOutput.close();

    return 0;
}

/*my_type**& myInverse(my_type**& A1, const int& size) {
    my_type m;
    my_type k;
    bool fl;
    bool noProblems = true;
    my_type dub;
    my_type* dub1 = new my_type;
    std::tuple<int, int> t;
    std::vector<std::tuple<int, int>> permutations = {};

    my_type** A = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        A[i] = new my_type[size];
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A[i][j] = A1[i][j];
        }
    }

    my_type** E = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        E[i] = new my_type[size]{};
    }

    for (int i = 0; i < size; i++) {
        E[i][i] = 1.0;
    }

    fl = true;
    for (int i = 0; i < size - 1; i++) {
        t = SearchMax(A, i, size);
        if (std::get<0>(t) != i) {
            std::swap(A[i], A[std::get<0>(t)]);
            std::swap(E[i], E[std::get<0>(t)]);
            /*
            dub1 = A[i];
            A[i] = A[std::get<0>(t)];
            A[i] = A[std::get<0>(t)];
            A[std::get<0>(t)] = dub1;


            dub1 = E[i];
            E[i] = E[std::get<0>(t)];
            E[std::get<0>(t)] = dub1;
        }

        if (EisDegenerate(A, i, size)) {
            fl = false;
        }

        if (fl == true) {
            m = A[i][i];
            for (int n = 0; n < size; n++) {
                E[i][n] /= m;
            }

            for (int j = i; j < size; j++) {
                A[i][j] /= m;
            }

            for (int i1 = i + 1; i1 < size; i1++) {

                k = A[i1][i];

                for (int n = 0; n < size; n++) {
                    E[i1][n] -= k * E[i][n];
                }

                for (int j1 = i; j1 < size; j1++) {
                    A[i1][j1] -= k * A[i][j1];
                }
            }
        }
        else {
            noProblems = false;
            break;
        }
    }

    t = SearchMax(A, size - 1, size);
    if (std::get<0>(t) != size - 1) {
        dub1 = A[size - 1];
        A[size - 1] = A[std::get<0>(t)];
        A[std::get<0>(t)] = dub1;
        dub1 = E[size - 1];
        E[size - 1] = E[std::get<0>(t)];
        E[std::get<0>(t)] = dub1;
    }
    if (EisDegenerate(A, size - 1, size)) {
        noProblems = false;
    }

    if (noProblems == true) {

        for (int n = 0; n < size; n++) {
            E[size - 1][n] /= A[size - 1][size - 1];
        }
        A[size - 1][size - 1] /= A[size - 1][size - 1];

        for (int i = size - 1; i > -1; i--) {
            for (int j = i - 1; j > -1; j--) {

                for (int n = 0; n < size; n++) {
                    E[j][n] -= A[j][i] * E[i][n];
                }
                A[j][i] = 0.0;
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (A[i][j] > epsilon && i != j) {
                    dub1 = A[i];
                    A[i] = A[std::get<0>(t)];
                    A[std::get<0>(t)] = dub1;
                    dub1 = E[i];
                    E[i] = E[std::get<0>(t)];
                    E[std::get<0>(t)] = dub1;
                }
            }
        }
    }
    permutations.clear();
    for (int i = 0; i < size; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    return E;
}*/

void ECalculations(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP, my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, my_type*& BF, const int& size) {
    
    my_type** T = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        T[i] = new my_type[size]{};
    }

    my_type** Q = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        Q[i] = new my_type[size]{};
    }

    my_type** R = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        R[i] = new my_type[size]{};
    }

    my_type** Buf = new my_type * [size];
    for (int i = 0; i < size; ++i) {
        Buf[i] = new my_type[size]{};
    }

    std::tuple<my_type, my_type, my_type> d11;
    std::tuple<my_type, my_type> c11;
    my_type cond;
    bool noProblem;
    my_type dx;
    my_type db;
    srand(3);

    for (int i = 0; i < size; i++) {
        dB[i] = B[i] * 0.01 - (my_type)(rand()) / RAND_MAX ;
    }

    EGaussMethod(A, B, X1, size, noProblem);

    if (noProblem) {

        EGaussMethod(A, dB, X2, size, noProblem);
        d11 = Ediscrepancy(A, B, X1,BF, size);

        EMatrixInverse(A,IA,T,Q,R,Buf,B,size);
        EMatrixMult(A, IA, size);

        c11 = EConditionNumber(A, IA, size);
        if (noProblem) {
            cond = EmarkCondition(X2, X1, dB, B, dx, db, size);
            for (int i = 0; i < 20; i++) {
                for (int i = 0; i < size; i++) {
                    dB[i] = B[i]*(my_type)(rand()) / RAND_MAX;
                }
                if (cond < EmarkCondition(X2, X1, dB, B, dx, db, size)) {
                    cond = EmarkCondition(X2, X1, dB, B, dx, db, size);
                }
            }
            EWriteData(fileNameA, fileNameIA, fileNameB,fileNameX1, fileNameX2, fileNameP,A, IA, X1, X2,
                B, std::get<0>(d11), std::get<1>(d11), std::get<2>(d11), cond, std::get<0>(c11),
                std::get<1>(c11), db, dx, size);
        }
       // for (int i = 0; i < size; i++) {
       //     dB[i] = B[i]/**0.001*/-  (my_type)(rand()) / RAND_MAX * 2;
       // }
    }
    else {
        std::cout << "In "+ fileNameA + " matrix A is degenerate " << std::endl;
    }
    for (int i = 0; i < size; ++i)
    {
        delete[] T[i];
        delete[] Q[i];
        delete[] R[i];
        delete[] Buf[i];
    }
    delete[] T;
    delete[] Q;
    delete[] R;
    delete[] Buf;
}


int EConditionNumberQR(my_type**& matrixR, my_type**& matrixT, my_type*& vector, const size_t column, const size_t n)
{
    //std::cout << "goes CONDNUM stage: " << column << std::endl;
    //WriteMatrix("matrixR: ", matrixR, n);
    size_t maxNumber = column;
    my_type maxValue = matrixR[column][column];
    for (int i = column; i < n; ++i)
    {
        //std::cout << "if " << fabs(matrixR[i][column]) << " > " << fabs(maxNumber) << std::endl;
        if (fabs(matrixR[i][column]) > fabs(maxValue))
        {
            maxValue = matrixR[i][column];
            maxNumber = i;
        }
    }

    if (maxNumber != column) //if diagonal element is not max
    {
        //std::cout << "stage " << column << " max number " << maxNumber << std::endl;
        std::swap(matrixR[column], matrixR[maxNumber]);
        std::swap(matrixT[column], matrixT[maxNumber]);
        std::swap(vector[column], vector[maxNumber]);
    }
    else
    {
        //std::cout << "maxValue:" << maxValue << std::endl;
        //std::cout << "maxNumber:" << maxNumber << std::endl;
    }

    return 0;
}

int EReverseMotion(my_type**& matrixR, my_type*& vectorX, my_type*& vectorB, const size_t n)
{
    // if(abs(matrixR[n-1][n-1]) <= epsilon)
    // {
    //     std::cout << matrixR[n-1][n-1] << std::endl;
    //     std::cout << "Matrix is singular\n";
    //     return 1;
    // }
    // else{
    vectorX[0] = 1;
    for (int i = n - 1; i >= 0; --i)
    {
        my_type sum = 0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += matrixR[i][j] * vectorX[j];
        }
        vectorX[i] = (vectorB[i] - sum) / matrixR[i][i];
    }
    // }
    return 0;
}

int EQRCalculations(my_type**& matrixA, my_type**& matrixT, my_type**& matrixQ, my_type**& matrixR, my_type*& vectorB,
    my_type*& vectorX, my_type**& matrixBuffer1, my_type**& matrixBuffer2,
    my_type*& vectorBStarred, const size_t n)
{
    EMatrixCopy(matrixR, matrixA, n);
    //QRDecomposer2(matrixA, matrixQ, matrixR, matrixBuffer1, matrixBuffer2, n);
    EQRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer1[0], matrixBuffer1[1], vectorB, n);
    EMatrixMult(matrixT, vectorB, vectorBStarred, n);
    EReverseMotion(matrixR, vectorX, vectorBStarred, n);
    return 0;
}

int EQRDecomposer(my_type**& matrixA, my_type**& matrixT, my_type**& matrixQ,
    my_type**& matrixR, my_type*& vectorBuffer1, my_type*& vectorBuffer2, my_type*& vectorB, const size_t n)
{
    EMatrixCopy(matrixR, matrixA, n);

    EGetMatrixI(matrixT, n);

    for (int i = 0; i < n - 1; ++i)
    {
        //тут выбор главного элемента
        EConditionNumberQR(matrixR, matrixT, vectorB, i, n);

        for (int j = i + 1; j < n; ++j)
        {
            my_type c = matrixR[i][i];
            my_type s = matrixR[j][i];

            my_type radical = 1 / sqrt(c * c + s * s);

            c *= radical;
            s *= radical;

            for (int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixR[i][k] + s * matrixR[j][k]; //matrixR[i][k]
                vectorBuffer2[k] = (-s) * matrixR[i][k] + c * matrixR[j][k]; //matrixR[j][k]
            }

            for (int k = 0; k < n; ++k)
            {
                matrixR[i][k] = vectorBuffer1[k];
                matrixR[j][k] = vectorBuffer2[k];
            }

            for (int k = 0; k < n; ++k)
            {
                vectorBuffer1[k] = c * matrixT[i][k] + s * matrixT[j][k];
                vectorBuffer2[k] = (-s) * matrixT[i][k] + c * matrixT[j][k];
            }

            for (int k = 0; k < n; ++k)
            {
                matrixT[i][k] = vectorBuffer1[k];
                matrixT[j][k] = vectorBuffer2[k];
            }
            //WriteMatrix("R", matrixR, n);
        }
    }

    EMatrixTranspose(matrixT, matrixQ, n);

    return 0;
}

int EMatrixInverseTR(my_type**& matrixT, my_type**& matrixR, my_type**& matrixInverted,
    my_type**& matrixBuffer, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixBuffer[0][j] = matrixT[j][i];
        }

        EReverseMotion(matrixR, matrixInverted[i], matrixBuffer[0], n);
    }

    EMatrixTranspose(matrixInverted, matrixBuffer, n);
    EMatrixCopy(matrixInverted, matrixBuffer, n);

    return 0;
}

int EMatrixInverse(my_type**& matrixA, my_type**& matrixInverted, my_type**& matrixT,
    my_type**& matrixQ, my_type**& matrixR, my_type**& matrixBuffer, my_type*& vectorB, const size_t n)
{
    EQRDecomposer(matrixA, matrixT, matrixQ, matrixR, matrixBuffer[0], matrixBuffer[1], vectorB, n);

    EMatrixInverseTR(matrixT, matrixR, matrixInverted, matrixBuffer, n);

    return 0;
}

int EMatrixTranspose(my_type**& matrixInit, my_type**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixResult[i][j] = matrixInit[j][i];
        }
    }

    return 0;
}

int EMatrixCopy(my_type**& matrixPaste, my_type**& matrixCopy, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matrixPaste[i][j] = matrixCopy[i][j];
        }
    }
    return 0;
}

int EGetMatrixI(my_type**& matrix, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                matrix[i][j] = 1;
            }
            else
            {
                matrix[i][j] = 0;
            }
        }
    }

    return 0;
}

int EMatrixMultV(my_type**& matrixA, my_type**& matrixB, my_type**& matrixResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            my_type sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum += matrixA[i][k] * matrixB[k][j];
            }
            matrixResult[i][j] = sum;
        }
    }

    return 0;
}//matrixA * matrixB

int EWriteMatrix(const std::string label, my_type**& matrix, const size_t n)
{
    std::cout << label << "\n";

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    return 0;
}

int EMatrixMult(my_type**& matrix, my_type*& vector, my_type*& vectorResult, const size_t n)
{
    for (int i = 0; i < n; ++i)
    {
        my_type sum = 0;
        for (int j = 0; j < n; ++j)
        {
            sum += matrix[i][j] * vector[j];
        }
        vectorResult[i] = sum;
    }

    return (0);
}//matrix * vector
