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

    // Initialize the solution based on the input file (or restart file)
    void initializeSolutionArrays() override;

    void solve() override;

    // void computeFluxes() override;
    // void spatialIntegration() override;

    // given a certain conservative solution, compute the dt corresponding to the CFL condition
    Matrix3D<FloatType> computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep);
    
    void updateMassFlows(const FlowSolution &solution);

    FlowSolution computeFluxResiduals(const FlowSolution solution, unsigned long int itCounter) const;

    void computeAdvectionResiduals(FluxDirection direction, const FlowSolution solution, unsigned long int itCounter, FlowSolution &residuals) const;

    void updateSolution(FlowSolution &solOld, FlowSolution &solNew, const FlowSolution &residuals, const FloatType &integrationCoeff, const Matrix3D<FloatType> &dt);

    void printInfoResiduals(FlowSolution &residuals, unsigned long int it) const;

    void printHeader() const;

    StateVector computeLogResidualNorm(const FlowSolution &residuals) const;

private:
    // // Add CEulerSolver-specific member variables here
    FlowSolution _conservativeVars;

    // std::vector<FloatType> _residualNorms;
    
    // // Example: time step
    // FloatType _dt;
};
