#include <iostream>
#include "CMesh.hpp"
#include "Config.hpp"
#include "CSolverBase.hpp"
#include "CSolverEuler.hpp"


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.ini>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    Config config(inputFile);
    
    CMesh mesh(config);    
    
    KindSolver kindSolver = config.getKindSolver();
    
    std::unique_ptr<CSolverBase> solver;
    
    if (kindSolver == KindSolver::EULER) {
       solver = std::make_unique<CSolverEuler>(config, mesh);
    }
    else {
        std::cerr << "Unsupported solver kind: " << static_cast<int>(kindSolver) << std::endl;
        return 1;
    }
    
    solver->solve();

    return 0;
}
