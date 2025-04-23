#include "../include/CEulerSolver.hpp"
#include <iostream> // Optional, for logging/debugging

CEulerSolver::CEulerSolver(Config& config, CMesh& mesh)
    : CSolverBase(config, mesh)  // Call base class constructor
{
    std::cout << "CEulerSolver initialized." << std::endl;
    instantiateSolutionArrays();
}

void CEulerSolver::instantiateSolutionArrays(){
    _conservativeVars.resize(_nPointsI, _nPointsJ, _nPointsK);
}