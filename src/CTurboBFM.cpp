#include <iostream>
#include "CMesh.hpp"
#include "Config.hpp"
#include "CEulerSolver.hpp"


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.ini>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    Config config(inputFile);
    CMesh mesh(config);    
    CEulerSolver solver(config, mesh);
    solver.solve();

    return 0;
}
