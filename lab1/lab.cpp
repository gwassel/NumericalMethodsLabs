#include <iostream>
#include <string>
#include "../nlohmann/json.hpp"

#define my_type double

using json = nlohmann::json;

int readInit(const std::string path);

int allocateMemory(my_type** &Matrix_A, my_type* &Matrix_b);

int readData(const std::string path);

int calculations(my_type** &Matrix_A, my_type* &Matrix_b);

int freeMemory(my_type** &Matrix_A, my_type* &Matrix_b);


int main()
{
    std::cout << "Hello, world!\n";
    return 0;
}
