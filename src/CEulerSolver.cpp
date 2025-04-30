#include "CEulerSolver.hpp"
#include "commonFunctions.hpp"
#include "types.hpp"
#include <iostream> // Optional, for logging/debugging
#include <iomanip>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

CEulerSolver::CEulerSolver(Config& config, CMesh& mesh)
    : CSolverBase(config, mesh)  // Call base class constructor
{
    initializeSolutionArrays();
}

void CEulerSolver::initializeSolutionArrays(){
    

    FloatType initMach = _config.getInitMachNumber();
    FloatType initTemperature = _config.getInitTemperature();
    FloatType initPressure = _config.getInitPressure();
    Vector3D initDirection = _config.getInitDirection();

    FloatType density {0.0}, totEnergy {0.0};
    Vector3D velocity {0.0, 0.0, 0.0};
    _fluid->computeInitFields(initMach, initTemperature, initPressure, initDirection, density, velocity, totEnergy);
    
    _conservativeVars.resize(_nPointsI, _nPointsJ, _nPointsK);
    for (int i=0; i<_nPointsI; i++) {
        for (int j=0; j<_nPointsJ; j++){
            for (int k=0; k<_nPointsK; k++){
                _conservativeVars._rho(i,j,k) = density;
                _conservativeVars._rhoU(i,j,k) = density * velocity.x();
                _conservativeVars._rhoV(i,j,k) = density * velocity.y();
                _conservativeVars._rhoW(i,j,k) = density * velocity.z();
                _conservativeVars._rhoE(i,j,k) = density * totEnergy;
            }
        }
    }

}

void CEulerSolver::solve(){
    size_t nIterMax = _config.getMaxIterations();
    Matrix3D<FloatType> timestep(_nPointsI, _nPointsJ, _nPointsK);                          // place holder for time step array
    std::vector<FloatType> timeIntegrationCoeffs = _config.getTimeIntegrationCoeffs();      // time integration coefficients (runge kutta)
    FlowSolution fluxResiduals(_nPointsI, _nPointsJ, _nPointsK);                            // place holder for flux residuals
    size_t updateMassFlowsFreq = _config.getMassFlowUpdateFrequency();                      // frequency to update the mass flows at the boundaries

    // time integration
    for (size_t it=0; it<nIterMax; it++){        
        FlowSolution solutionOld = _conservativeVars;                                       // place holder for the solution at the previous timestep
        if (it%updateMassFlowsFreq == 0) updateMassFlows(solutionOld);                      // compute the mass flows
        computeTimestepArray(solutionOld, timestep);                                        // compute the physical time step
        FlowSolution tmpSol = solutionOld;                                                  // place holder for the solution at the runge-kutta step
        
        // runge-kutta steps
        for (const auto &integrationCoeff: timeIntegrationCoeffs){
            fluxResiduals = computeFluxResiduals(tmpSol, it);
            updateSolution(solutionOld, tmpSol, fluxResiduals, integrationCoeff, timestep);
        }
        
        // update the solution and print information
        _conservativeVars = tmpSol;
        printInfoResiduals(fluxResiduals, it);
        printInfoMassFlows(it);
        if (it%updateMassFlowsFreq == 0) saveSolution(); 
    }
}

void CEulerSolver::printInfoResiduals(FlowSolution &residuals, size_t it) const {
    if (it == 0) {printHeader();}
    auto logRes = computeLogResidualNorm(residuals);
    printLogResiduals(logRes, it);
}

void CEulerSolver::printInfoMassFlows(size_t it) const {
    if (it%100 == 0){
        std::cout << "\nMASS FLOWS CHECK [kg/s]:\n";
        std::cout << "I_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::I_START) << std::endl;
        std::cout << "I_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::I_END) << std::endl;
        std::cout << "J_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::J_START) << std::endl;
        std::cout << "J_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::J_END) << std::endl;
        std::cout << "K_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::K_START) << std::endl;
        std::cout << "K_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::K_END) << std::endl << std::endl;
    }
}

void CEulerSolver::printLogResiduals(const StateVector &logRes, unsigned long int it) const {
    int col_width = 14;
    std::cout << std::fixed << std::setprecision(6); // Set fixed format with 6 decimals
    std::cout << "|" << std::setw(col_width) << std::setfill(' ') << std::left << it+1 << "|"
              << std::setw(col_width) << std::left << 0.0 << "|"
              << std::setw(col_width) << std::right << logRes[0] << "|"
              << std::setw(col_width) << std::right << logRes[1] << "|"
              << std::setw(col_width) << std::right << logRes[2] << "|"
              << std::setw(col_width) << std::right << logRes[3] << "|"
              << std::setw(col_width) << std::right << logRes[4] << "|"
              << std::endl;
}



StateVector CEulerSolver::computeLogResidualNorm(const FlowSolution &residuals) const {
    StateVector logResidualNorm{};
    for (int i=0; i<5; i++){
        auto residualNorm = residuals.norm(i);
        if (residualNorm >= 1E-12) {
            logResidualNorm[i] = std::log10(residuals.norm(i) / (_nPointsI*_nPointsJ*_nPointsK));
        } else{
            logResidualNorm[i] = 0.0;
        }
    }
    return logResidualNorm;
}

void CEulerSolver::printHeader() const {
    int col_width = 14;
    
    // Print the first separator (top border)
    std::cout << "|" << std::setw(col_width * 7 + 6) << std::setfill('-') << "" << "|" << std::endl;

    // Print the column headers
    std::cout << "|"
                << std::setw(col_width) << std::setfill(' ') << std::left << "Iteration" << "|"
                << std::setw(col_width) << std::left << "Time[s]" << "|"
                << std::setw(col_width) << std::right << "rms[Rho]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoU]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoV]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoW]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoE]" << "|"
                << std::endl;

    // Print the last separator (bottom border)
    std::cout << "|" << std::setw(col_width * 7 + 6) << std::setfill('-') << "" << "|" << std::endl;
}

Matrix3D<FloatType> CEulerSolver::computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep){
    FloatType cflMax = _config.getCFL();
    Vector3D iEdge, jEdge, kEdge;
    Vector3D iDir, jDir, kDir;
    Vector3D velocity;
    StateVector primitive, conservative;
    FloatType soundSpeed;
    std::array<FloatType, 3> dtEdge;
    FloatType dtMin;

    for (size_t i=0; i<_nPointsI; i++) {
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                _mesh.getElementEdges(i, j, k, iEdge, jEdge, kEdge);
                iDir = iEdge / iEdge.magnitude();
                jDir = jEdge / jEdge.magnitude();
                kDir = kEdge / kEdge.magnitude();
                conservative = solution.at(i,j,k);
                primitive = getEulerPrimitiveFromConservative(conservative);
                velocity(0) = primitive[1];
                velocity(1) = primitive[2];
                velocity(2) = primitive[3];
                soundSpeed = _fluid->computeSoundSpeed_rho_u_et(primitive[0], velocity, primitive[4]);

                dtEdge[0] = iEdge.magnitude() / (std::abs(velocity.dot(iDir))+soundSpeed);
                dtEdge[1] = jEdge.magnitude() / (std::abs(velocity.dot(jDir))+soundSpeed);
                dtEdge[2] = kEdge.magnitude() / (std::abs(velocity.dot(kDir))+soundSpeed);

                dtMin = *std::min_element(dtEdge.begin(), dtEdge.end());
                timestep(i,j,k) = dtMin*cflMax;
            }
        }
    }
    return timestep;
}


void CEulerSolver::updateMassFlows(const FlowSolution&solution){
    std::array<BoundaryIndices, 6> bcIndices {BoundaryIndices::I_START,
                                              BoundaryIndices::I_END,
                                              BoundaryIndices::J_START,
                                              BoundaryIndices::J_END,
                                              BoundaryIndices::K_START,
                                              BoundaryIndices::K_END};
    
    for (auto& bcIndex: bcIndices){
        Matrix2D<Vector3D> surface = _mesh.getMeshBoundary(bcIndex);
        Matrix2D<FloatType> rhoUX = (_conservativeVars._rhoU).getBoundarySlice(bcIndex);
        Matrix2D<FloatType> rhoUV = (_conservativeVars._rhoV).getBoundarySlice(bcIndex);
        Matrix2D<FloatType> rhoUW = (_conservativeVars._rhoW).getBoundarySlice(bcIndex);
        _massFlows[bcIndex] = computeSurfaceIntegral(surface, rhoUX, rhoUV, rhoUW);
    }
    
    
}

FlowSolution CEulerSolver::computeFluxResiduals(const FlowSolution& solution, size_t it) const {
    FlowSolution residuals(_nPointsI, _nPointsJ, _nPointsK);
    computeAdvectionResiduals(FluxDirection::I, solution, it, residuals);
    computeAdvectionResiduals(FluxDirection::J, solution, it, residuals);
    if (_topology==Topology::THREE_DIMENSIONAL || _topology==Topology::AXISYMMETRIC){
        computeAdvectionResiduals(FluxDirection::K, solution, it, residuals);
    }
    return residuals;
}

void CEulerSolver::computeAdvectionResiduals(FluxDirection direction, const FlowSolution& solution, size_t itCounter, FlowSolution &residuals) const {
    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);
    const Matrix3D<Vector3D>& midPoints = _mesh.getMidPoints(direction);
    
    BoundaryIndices boundaryStart, boundaryEnd;
    if (direction==FluxDirection::I){
        boundaryStart = BoundaryIndices::I_START;
        boundaryEnd = BoundaryIndices::I_END;
    }
    else if (direction==FluxDirection::J){
        boundaryStart = BoundaryIndices::J_START;
        boundaryEnd = BoundaryIndices::J_END;
    }
    else {
        boundaryStart = BoundaryIndices::K_START;
        boundaryEnd = BoundaryIndices::K_END;
    }
    
    StateVector Uinternal{}, Uleft{}, Uright{}, Uleftleft {}, Urightright {}, flux {};
    Vector3D surface {}, midPoint {};

    size_t ni = surfaces.sizeI(); 
    size_t nj = surfaces.sizeJ(); 
    size_t nk = surfaces.sizeK();

    for (size_t iFace = 0; iFace < ni; ++iFace) {
        for (size_t jFace = 0; jFace < nj; ++jFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
                size_t dirFace = 0;
                size_t stopFace = 0;
                switch (direction)
                {
                case (FluxDirection::I):
                    dirFace = iFace;
                    stopFace = ni-1;
                    break;
                case (FluxDirection::J):
                    dirFace = jFace;
                    stopFace = nj-1;
                    break;
                case (FluxDirection::K):
                    dirFace = kFace;
                    stopFace = nk-1;
                    break;
                default:
                    throw std::runtime_error("Invalid FluxDirection.");
                }
                
                // fluxes calculation here, also boundary conditions.
                if (dirFace == 0) {
                    Uinternal = solution.at(iFace, jFace, kFace);
                    surface = -surfaces(iFace, jFace, kFace);
                    midPoint = midPoints(iFace, jFace, kFace);
                    flux = _boundaryConditions.at(boundaryStart)->computeBoundaryFlux(Uinternal, surface, midPoint);
                    residuals.add(iFace, jFace, kFace, flux * surface.magnitude());
                } else if (dirFace == stopFace) {
                    Uinternal = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
                    surface = surfaces(iFace, jFace, kFace);
                    midPoint = midPoints(iFace, jFace, kFace);
                    flux = _boundaryConditions.at(boundaryEnd)->computeBoundaryFlux(Uinternal, surface, midPoint);
                    residuals.add(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2], flux * surface.magnitude());
                } else {
                    // internal flux calculation
                    Uleft = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
                    Uright = solution.at(iFace, jFace, kFace);

                    if (dirFace==1){
                        Uleftleft = Uleft;
                        Urightright = solution.at(iFace+1*stepMask[0], jFace+1*stepMask[1], kFace+1*stepMask[2]);
                    }
                    else if (dirFace==stopFace-1){
                        Uleftleft = solution.at(iFace-2*stepMask[0], jFace-2*stepMask[1], kFace-2*stepMask[2]);
                        Urightright = Uright;
                    }
                    else {
                        Uleftleft = solution.at(iFace-2*stepMask[0], jFace-2*stepMask[1], kFace-2*stepMask[2]);
                        Urightright = solution.at(iFace+1*stepMask[0], jFace+1*stepMask[1], kFace+1*stepMask[2]);
                    }

                    surface = surfaces(iFace, jFace, kFace);
                    flux = _advectionScheme->computeFlux(Uleftleft, Uleft, Uright, Urightright, surface);
                    residuals.add(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2], flux * surface.magnitude());
                    residuals.subtract(iFace, jFace, kFace, flux * surface.magnitude());
                }

            
            }
        }
    }
    
}

void CEulerSolver::updateSolution(const FlowSolution &solOld, FlowSolution &solNew, const FlowSolution &residuals, const FloatType &integrationCoeff, const Matrix3D<FloatType> &dt){
    for (size_t i = 0; i < _nPointsI; i++) {
        for (size_t j = 0; j < _nPointsJ; j++) {
            for (size_t k = 0; k < _nPointsK; k++) {
                solNew.set(i, j, k, solOld.at(i, j, k) - residuals.at(i, j, k) * integrationCoeff * dt(i, j, k) / _mesh.getVolume(i,j,k));
            }
        }
    }
}


void CEulerSolver::saveSolution() const {
    

}