#include <iostream>
#include "../nlohmann/json.hpp"


using json = nlohmann::json;

int readInit();

int allocateMemory();

int readData();

int calculations();

int freeMemory();


int main()
{
    std::cout << "Hello, world!\n";
    return 0;
}
