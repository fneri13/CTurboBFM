#include "../../include/CMesh.hpp"
#include "../../include/Config.hpp"
#include "../../include/types.hpp"
#include "../../include/CEulerSolver.hpp"
#include <iostream>

int main() {
    Config config("input.ini");
    CMesh mesh(config);    
    // KindSolver solverType = config.getKindSolver();
    CEulerSolver solver(config, mesh);

    return 0;
}
