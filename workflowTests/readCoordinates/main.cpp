#include "../../include/CMesh.hpp"
#include "../../include/Config.hpp"
#include <iostream>

int main() {
    std::cout << "Reading Config file..." << std::endl;
    Config config("input.ini");
    CMesh mesh(config);
    return 0;
}
