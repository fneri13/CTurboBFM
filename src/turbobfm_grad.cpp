#include <iostream>
#include "mesh.hpp"
#include "config.hpp"
#include "solver_base.hpp"
#include "solver_euler.hpp"


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.ini>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    Config config(inputFile);
    
    Mesh mesh(config);    
    
    KindSolver kindSolver = config.getKindSolver();
    
    std::unique_ptr<SolverBase> solver;
    
    if (kindSolver == KindSolver::EULER) {
       solver = std::make_unique<SolverEuler>(config, mesh);
    }
    else {
        std::cerr << "Unsupported solver kind: " << static_cast<int>(kindSolver) << std::endl;
        return 1;
    }
    
    solver->writeSolution(0, true);

    return 0;
}
