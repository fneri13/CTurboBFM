#include "../../include/CMesh.hpp"
#include <iostream>

int main() {
    std::cout << "Starting mesh generation..." << std::endl;
    CMesh mesh("coordinates.csv");
    return 0;
}
