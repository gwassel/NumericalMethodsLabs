#include "funcs.hpp"

const double epsilon = 1e-12;

void printEquation(my_type**& A, my_type*& B, const int& size, std::string s) {
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

void printMatrix(my_type**& A, const int& size, std::string s) {
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void printVector(my_type*& B, const int& size, std::string s){
    std::cout << s << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << B[i] << std::endl;
    }
    std::cout << std::endl;
}

void GaussMethod(my_type** A1, my_type* B, my_type*& X, const int& size, bool& flag) {
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
        if (MatrixIsPrepared(A, B, i, permutations, size)) {
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
    if (MatrixIsPrepared(A,B,size - 1, permutations, size) == false) {
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
        DiagonalizeEquation(A,B,X,size,permutations);
        for (int i = 0; i < size; i++) {
            dub = X[i];
            X[i] = B[i];
            B[i] = dub;
        }
    }
    else {
        //printMatrix(A,size,"Deg");
    }
    flag = noProblems;
    permutations.clear();
}

bool MatrixIsPrepared(my_type**& A, my_type*& B, const int& i, std::vector<std::tuple<int, int>>& permutations, const int& size) {
    std::tuple<int, int> t;
    t = SearchMax(A,i,size);

    if (std::get<0>(t) != i) {
        SwapRows(A,B,i, std::get<0>(t));
    }
    if (std::get<1>(t) != i) {
        SwapColomns(A,i, std::get<1>(t),size);
        permutations.push_back(std::make_tuple(i, std::get<1>(t)));
    }
    if (isDegenerate(A,i,size) == false) {
        return true;
    }
    else return false;
}

void DiagonalizeEquation(my_type**& A, my_type*& B, my_type*& X, const int& size, std::vector<std::tuple<int, int>>& permutations) {
    for (int i = permutations.size() - 1; i > -1; i--) {
        SwapColomns(A, std::get<0>(permutations[i]), std::get<1>(permutations[i]), size);
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (A[i][j] > epsilon && i != j) { //äðóãîå ñðàâíåíèå ñ íóëåì
                SwapRows(A, B, i, j);
            }
        }
    }
    permutations.clear();
}

std::tuple<int, int> SearchMax(my_type**& A, const int& currentMinor, const int& size) {
    int max_i = currentMinor;
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
    return std::make_tuple(max_i, max_j);
}

void SwapRows(my_type**& A, my_type*& B, const int& i1, const int& i2) {
    my_type* dub1 = new my_type;
    dub1 = A[i1];
    A[i1] = A[i2];
    A[i2] = dub1;

    my_type dub2 = B[i1];
    B[i1] = B[i2];
    B[i2] = dub2;
    //std::cout << "SwapRows " << i1 << " " << i2 << std::endl;
    //printEquation("Swap Rows!");
}

void SwapColomns(my_type**& A, const int& j1, const int& j2, const int& size) {
    my_type dub;
    for (int i = 0; i < size; i++) {
        dub = A[i][j1];
        A[i][j1] = A[i][j2];
        A[i][j2] = dub;
    }

    //std::cout << "SwapColomns "<<j1<<" "<<j2 << std::endl;
    //printEquation("Swap Colomns!");
}

bool isDegenerate(my_type**& A, const int& i,const int& size) {
   /* int k = 0;
    for (int j = 0; j < size; j++) {
        if (fabs(A[i][j] / A[i][j]) != fabs(A[i][j] / A[i][j])) {
            k++;
        }
    }*/
    if (/*k == size*//*fabs(A[i][i] / A[i][i]) != fabs(A[i][i] / A[i][i])*/fabs(A[i][i])<epsilon) {
        return true;
    }
    return false;
}


my_type CubicVectorNorm(my_type* p, const int& size) {
    my_type sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(p[i]);
    }
    return sum;
}

my_type CubicMatrixNorm(my_type** p, const int& size) {
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

my_type OctahedralVectorNorm(my_type* p, const int& size) {
    my_type max = fabs(p[0]);
    for (int i = 1; i < size; i++) {
        if (max < fabs(p[i]))
            max = fabs(p[i]);
    }
    return max;
}

my_type OctahedralMatrixNorm(my_type** p, const int& size) {
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

std::tuple<my_type, my_type> ConditionNumber(my_type**& A, my_type**& IA, const int& size) {

    my_type c;

    c = CubicMatrixNorm(IA, size) * CubicMatrixNorm(A, size);
    //std::cout << "Condition number with cubic matrix norm:" << c << std::endl;

    c = OctahedralMatrixNorm(IA, size) * OctahedralMatrixNorm(A, size);
    //std::cout << "Condition number with octahedral matrix norm:" << c << std::endl;
    //std::cout << std::endl;
    return std::make_tuple(CubicMatrixNorm(IA, size) * CubicMatrixNorm(A, size), OctahedralMatrixNorm(IA, size) * OctahedralMatrixNorm(A, size));
}

std::tuple<my_type,my_type,my_type> discrepancy(my_type**& A, my_type*& B, my_type*& X, const int& size) {
    my_type* D = new my_type[size];

    for (int i = 0; i < size; i++)
    {
        D[i] = 0;
        for (int j = 0; j < size; j++)
        {
            D[i] += A[i][j] * X[j];
        }
    }

    my_type sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(D[i] - B[i], 2);
    }

    //std::cout << "Discrepancy(euclid) " << sqrt(sum) << std::endl;
   // std::cout << "Discrepancy(cubic) " << CubicVectorNorm(D, size) << std::endl;
    //std::cout << "Discrepancy(octahedral) " << OctahedralVectorNorm(D, size) << std::endl;
   // std::cout << std::endl;
    return std::make_tuple(sqrt(sum), CubicVectorNorm(D, size), OctahedralVectorNorm(D, size));
}

int FreeMemory(my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, const int& n)
{
    delete[] B;
    delete[] dB;
    delete[] X1;
    delete[] X2;

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
        delete[] IA[i];
    }
    delete[] A; 
    delete[] IA;

    return 0;
}

my_type*& minusVectors(my_type*& X1, my_type*& X2, const int& size) {
    for (int i = 0; i < size; i++) {
        X1[i] -= X2[i];
    }
    return X1;
}

my_type markCondition(my_type*& X1, my_type*& X2, my_type*& B1, my_type*& B2, my_type& dx, my_type& db, const int& size) {
    
    my_type normX1 = CubicVectorNorm(X1,size);
    my_type normX2 = CubicVectorNorm(X2, size);
    my_type normB1 = CubicVectorNorm(B1,size);
    my_type normB2 = CubicVectorNorm(B2,size);
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


int WriteData(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP,my_type**& matrixA, my_type**& matrixIA,
    my_type*& vectorX1, my_type*& vectorX2, my_type*& vectorB, const my_type dis1,
    const my_type dis2, const my_type dis3, const my_type condb, const my_type condNum1,
    const my_type condNum2, const my_type db, const my_type dx, const int& n)
{
    WriteMatrix(fileNameA, "matrix A", matrixA, n);
    WriteMatrix(fileNameIA, "matrix IA", matrixIA, n);

    WriteVector(fileNameB, "vector B", vectorB, n);
    WriteVector(fileNameX1, "vector X1", vectorX1, n);
    WriteVector(fileNameX2, "vector X2", vectorX2, n);

    WriteParameters(fileNameP, dis1, dis2, dis3, condb, condNum1, condNum2, db, dx);
    return 0;
}

int WriteParameters(std::string fileNameOutput, const my_type dis1, 
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

void GaussCalculations(std::string fileNameA, std::string fileNameIA, std::string fileNameB, std::string fileNameX1,
    std::string fileNameX2, std::string fileNameP, my_type**& A, my_type**& IA, my_type*& B, my_type*& dB, my_type*& X1, my_type*& X2, const int& size) {
    
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

    GaussMethod(A, B, X1, size, noProblem);

    if (noProblem) {

        GaussMethod(A, dB, X2, size, noProblem);
        d11 = discrepancy(A, B, X1, size);
        //IA = InverseMatrix(A, size);
        //c11 = ConditionNumber(A, IA, size);
        if (noProblem) {
            cond = markCondition(X2, X1, dB, B, dx, db, size);
            for (int i = 0; i < 20; i++) {
                for (int i = 0; i < size; i++) {
                    dB[i] = B[i]*(my_type)(rand()) / RAND_MAX;
                }
                if (cond < markCondition(X2, X1, dB, B, dx, db, size)) {
                    cond = markCondition(X2, X1, dB, B, dx, db, size);
                }
            }
            WriteData(fileNameA, fileNameIA, fileNameB,fileNameX1, fileNameX2, fileNameP,A, IA, X1, X2,
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
}