#pragma once
#include "types.hpp"
#include "CSolverBase.hpp"
#include "COutputBase.hpp"
#include "COutputCSV.hpp"
#include "CSourceBFMBase.hpp"
#include "CSourceBFMHall.hpp"
#include "CSourceBFMThollet.hpp"
#include "CSourceBFMLiftDrag.hpp"
#include "CSourceBFMFrozenForce.hpp"
#include "CSourceBFMFrozenGradient.hpp"
#include "CSourceBFMNeri.hpp"
#include "CSourceBFMChima.hpp"
#include "CSourceBFMLamprakis.hpp"
#include <unordered_map>

#ifdef _OPENMP 
    #include <omp.h>
#endif

/**
 * \brief     Solver for the Euler equations.
 * \details   Implements the specific logic for solving the compressible Euler equations using a finite volume method.
 * \author    Francesco Neri
 */
class CSolverEuler : public CSolverBase {
public:

    /**
     * @brief Constructs the object with given references.
     * @param config The configuration object.
     * @param mesh The mesh object.
     */
    CSolverEuler(Config& config, CMesh& mesh);

    ~CSolverEuler() override = default;

    /** @brief Solves the Euler equations with explicit time integration.*/
    void solve() override;

    virtual void checkConvergence(bool &exitLoop, bool &isSteady) const override;

    void writeSolution(size_t iterationCounter, bool alsoGradients=false) override;
    

private:
    
    FlowSolution _conservativeVars; // conservative variables solution

    std::map<SolutionNames, Matrix3D<Vector3D>> _solutionGrad;
    
    std::unique_ptr<COutputBase> _output;
    
    std::unique_ptr<CSourceBFMBase> _bfmSource;
    
    std::vector<StateVector> _logResiduals;
    
    Matrix3D<Vector3D> _inviscidForce, _viscousForce;
    
    Matrix3D<FloatType> _deviationAngle;

    std::vector<unsigned int> _monitorPoints_idxI, _monitorPoints_idxJ, _monitorPoints_idxK;
    unsigned int _numberMonitorPoints = 0;

    /** @brief Initialize the solution arrays.*/
    void initializeSolutionArrays() override;


    /** \brief compute the time-step array corresponding to the CFL condition 
     * \param[in] solution The conservative variables solution.
     * \param[out] timestep The time-step array reference.
    */
    void computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep);


    /** @brief update the mass flows at the boundaries 
     * @param[in] solution The conservative variables solution
    */
    void updateMassFlows(const FlowSolution &solution);


    /** @brief update the turbo performance data 
     * @param[in] solution The conservative variables solution
    */
    void updateTurboPerformance(const FlowSolution &solution);

    /** @brief update the monitor points data 
     * @param[in] solution The conservative variables solution
    */
    void updateMonitorPoints(const FlowSolution &solution);

    /** @brief initialize the monitor points data 
    */
    void initializeMonitorPoints();


    /** @brief compute the global residual, defined as V*du/dt + R = 0 
     * @param[in] solution The conservative variables solution
     * @param[in] itCounter The iteration counter
     * @return The global residual 3D matrix
    */
    void computeResiduals(const FlowSolution& solution, const std::map<SolutionNames, Matrix3D<Vector3D>>& solutionGrad, 
                            const size_t itCounter, const FloatType timePhysical, FlowSolution &residuals);


    /** @brief compute the residuals contribution from advection fluxes
     * @param[in] direction The flux direction
     * @param[in] solution The conservative variables solution
     * @param[in] itCounter The iteration counter
     * @param[out] residuals The residuals 3D matrix
    */
    void computeAdvectionFlux(FluxDirection direction, const FlowSolution& solution, size_t itCounter, FlowSolution &residuals) const;


    /** @brief compute the residuals contribution from viscous fluxes
     * @param[in] direction The flux direction
     * @param[in] solution The conservative variables solution
     * @param[in] itCounter The iteration counter
     * @param[out] residuals The residuals 3D matrix
    */
    void computeViscousFlux(FluxDirection direction, const FlowSolution& solution, const std::map<SolutionNames, Matrix3D<Vector3D>>& solutionGrad , size_t itCounter, FlowSolution &residuals) const;

    /** @brief update the conservative variables solution
     * @param[in] solOld The old conservative variables solution
     * @param[out] solNew The new conservative variables solution
     * @param[in] residuals The residuals 3D matrix
     * @param[in] integrationCoeff The integration coefficient (runge-kutta)
     * @param[in] dt The time step array
    */
    void updateSolution(const FlowSolution &solOld, FlowSolution &solNew, const FlowSolution &residuals, const FloatType &integrationCoeff, const Matrix3D<FloatType> &dt);

    /** @brief post process the solution, to account for periodic boundaries
     * @param[out] sol The new conservative variables solution
    */
    void postprocessSolution(FlowSolution &sol);

    /** @brief update the gradient of the solution
    */
    void computeSolutionGradient(FlowSolution &sol, std::map<SolutionNames, Matrix3D<Vector3D>> &solutionGrad);


    /** @brief print information about the residuals
     * @param[in] residuals The global residuals 3D matrix
     * @param[in] it The iteration counter
    */
    void printInfoResiduals(FlowSolution &residuals, size_t it);


    /** @brief initialize the solution from scratch*/
    void initializeSolutionFromScratch();


    /** @brief initialize the solution from restart file*/
    void initializeSolutionFromRestart();


    /** @brief print information about the mass flows
     * @param[in] it The iteration counter
    */
    void printInfoMassFlows(size_t it) const;


    /** @brief print information about the turbo performance
     * @param[in] it The iteration counter
    */
    void printInfoTurboPerformance(size_t it) const;


    /** @brief print the header of the residuals table*/
    void printHeader() const;


    /** @brief compute the residuals log10
     * @param[in] residuals The global residuals 3D matrix
     * @return The log10 of the normalized norm of the residuals
    */
    StateVector computeLogResidualNorm(const FlowSolution &residuals) const;


    /** @brief print the residuals log10*/
    void printLogResiduals(const StateVector &logRes, unsigned long int it) const;


    /** @brief write the residuals log10 to a csv file*/
    void writeLogResidualsToCSV() const;


    /** @brief write the turbo performance to a csv file*/
    void writeTurboPerformanceToCSV() const;

    /** @brief write the turbo performance to a csv file*/
    void writeMonitorPointsToCSV() const;


    /** @brief update the radial profiles (pressure) at the outlet, needed by boundary conditions.
     * @param[in] solution The current conservative variables solution
    */
    void updateRadialProfiles(FlowSolution &solution);


    /**
     * \brief Compute the residuals contribution from the source terms. Store the BFM inviscid and viscous forces at the same time
     * \param[in] solution The conservative variables solution
     * \param[in] itCounter The iteration counter
     * \param[out] residuals The reference to the residuals matrix
     * \param[out] inviscidForce The reference to the inviscid force matrix to update
     * \param[out] viscousForce The reference to the viscous force matrix to update
     */
    void computeSourceTerms(const FlowSolution& solution, const std::map<SolutionNames, Matrix3D<Vector3D>>& solutionGrad, 
                            size_t itCounter, FlowSolution &residuals, Matrix3D<Vector3D> &inviscidForce, 
                            Matrix3D<Vector3D> &viscousForce, Matrix3D<FloatType> &deviationAngle, FloatType timePhysical);


    /** @brief read the restart file and store data in the reference arguments
     * @param[in] restartFileName The restart file name
     * @param[out] NI The number of points in the i direction
     * @param[out] NJ The number of points in the j direction
     * @param[out] NK The number of points in the k direction
     * @param[out] inputDensity The density matrix
     * @param[out] inputVelX The x velocity matrix
     * @param[out] inputVelY The y velocity matrix
     * @param[out] inputVelZ The z velocity matrix
     * @param[out] inputTemperature The temperature matrix
    */
    void readRestartFile(const std::string &restartFileName, size_t &NI, size_t &NJ, size_t &NK,
                        Matrix3D<FloatType> &inputDensity, Matrix3D<FloatType> &inputVelX, Matrix3D<FloatType> &inputVelY, Matrix3D<FloatType> &inputVelZ, 
                        Matrix3D<FloatType> &inputTemperature, Matrix3D<Vector3D> &inputForceViscous, Matrix3D<Vector3D> &inputForceInviscid);
    

    /** @brief initialize the solution from restart file in axysymmetric mode (NI and NJ must be equal to input file)
     * @param[out] inputDensity The density matrix
     * @param[out] inputVelX The x velocity matrix
     * @param[out] inputVelY The y velocity matrix
     * @param[out] inputVelZ The z velocity matrix
     * @param[out] inputTemperature The temperature matrix
    */
    void axisymmetricRestart(Matrix3D<FloatType> &inputDensity, Matrix3D<FloatType> &inputVelX, Matrix3D<FloatType> &inputVelY, Matrix3D<FloatType> &inputVelZ, Matrix3D<FloatType> &inputTemperature);


    /** @brief initialize the solution from restart file in standard mode (NI,NJ,NK must be equal to input file)
     * @param[out] inputDensity The density matrix
     * @param[out] inputVelX The x velocity matrix
     * @param[out] inputVelY The y velocity matrix
     * @param[out] inputVelZ The z velocity matrix
     * @param[out] inputTemperature The temperature matrix
    */
    void standardRestart(Matrix3D<FloatType> &inputDensity, Matrix3D<FloatType> &inputVelX, Matrix3D<FloatType> &inputVelY, 
                         Matrix3D<FloatType> &inputVelZ, Matrix3D<FloatType> &inputTemperature, Matrix3D<Vector3D> &inputForceViscous, 
                         Matrix3D<Vector3D> &inputForceInviscid);
    

    void computeGradient(const Matrix3D<FloatType> &var, Matrix3D<Vector3D> &grad) const;

    StateVector computeGongSource(const FloatType& radius, const FloatType& theta, const FloatType& omega, const StateVector& primitive, 
                            const Vector3D& densityGrad, const Vector3D& velXGrad, const Vector3D& velYGrad, const Vector3D& velZGrad, 
                            const Vector3D& totEnergyGrad, const FloatType& volume) const;
    
    void setMomentumSolution(FlowSolution &residuals, const BoundaryIndices &boundaryIndices, const Vector3D &wallVelocity) const;
    
    /** @brief update the radial profiles (pressure) at the outlet, needed by boundary conditions
     * and enforce strong conditions for no-slip walls
     * @param[in,out] solution The reference to the conservative solution matrix
    */
    void preprocessSolution(FlowSolution &solution);


    StateVector computeViscousFlux(const StateVector& conservative, const Vector3D& velXGrad, const Vector3D& velYGrad, 
                                   const Vector3D& velZGrad, const Vector3D& tempGrad, const Vector3D& surface) const;

};
