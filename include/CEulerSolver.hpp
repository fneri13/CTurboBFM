#pragma once
#include "types.hpp"
#include "CSolverBase.hpp"

/**
 * \brief     Solver for the Euler equations.
 * \details   Implements the specific logic for solving the compressible Euler equations
 *            using a finite volume method.
 * \author    Francesco Neri
 */
class CEulerSolver : public CSolverBase {
public:
    // Constructor
    CEulerSolver(Config& config, CMesh& mesh);

    // Destructor
    ~CEulerSolver() override = default;

    // void initializeSolutionArrays() override;
    // void solve() override;
    // void computeFluxes() override;
    // void spatialIntegration() override;
    // void computeTimestepArray() override;
    // void checkConvergence() override;
    // void printInfoResiduals() override;
    // void saveSolution() override;
    // void computeGradient() override;
    // void interpolateScalar() override;

private:
    // // Add CEulerSolver-specific member variables here
    EulerConservativeVariables _conservativeVars;

    // std::vector<FloatType> _residualNorms;
    
    // // Example: time step
    // FloatType _dt;
};
