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
    CEulerSolver(Config& config, CMesh& mesh);

    ~CEulerSolver() override = default;

    // Initialize the solution based on the input file (or restart file)
    void initializeSolutionArrays() override;

    // start the explicit time-marching solution
    void solve() override;

    // given a certain conservative solution, compute the dt array corresponding to the CFL condition
    Matrix3D<FloatType> computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep);
    
    void updateMassFlows(const FlowSolution &solution);

    // compute the residuals of all fluxes
    FlowSolution computeFluxResiduals(const FlowSolution solution, unsigned long int itCounter) const;

    // compute the residuals for the advection fluxes
    void computeAdvectionResiduals(FluxDirection direction, const FlowSolution solution, unsigned long int itCounter, FlowSolution &residuals) const;

    // update the solution new with the onl and residuals
    void updateSolution(FlowSolution &solOld, FlowSolution &solNew, const FlowSolution &residuals, const FloatType &integrationCoeff, const Matrix3D<FloatType> &dt);

    // print information about the residuals
    void printInfoResiduals(FlowSolution &residuals, unsigned long int it) const;

    // print header at first iteration
    void printHeader() const;

    // compute the log10 of the norm of residuals
    StateVector computeLogResidualNorm(const FlowSolution &residuals) const;

    // print the residuals log10
    void printLogResiduals(const StateVector &logRes, unsigned long int it) const;

private:
    
    FlowSolution _conservativeVars; // conservative variables solution

};
