#include "../include/CEulerSolver.hpp"
#include "../include/commonFunctions.hpp"
#include "../include/types.hpp"
#include <iostream> // Optional, for logging/debugging
#include <iomanip>
#include <algorithm>
#include <vector>

CEulerSolver::CEulerSolver(Config& config, CMesh& mesh)
    : CSolverBase(config, mesh)  // Call base class constructor
{
    std::cout << "CEulerSolver initialized." << std::endl;

    _conservativeVars.resize(_nPointsI, _nPointsJ, _nPointsK);

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
    unsigned long int nIterMax = _config.getMaxIterations();
    // FloatType physicalTime {0.0};
    Matrix3D<FloatType> dt(_nPointsI, _nPointsJ, _nPointsK);
    std::vector<FloatType> timeIntegrationCoeffs = _config.getTimeIntegrationCoeffs();
    FlowSolution fluxResiduals(_nPointsI, _nPointsJ, _nPointsK);
    FlowSolution solutionOld = _conservativeVars;
    FlowSolution solutionNew = _conservativeVars;

    for (unsigned long int it=0; it<nIterMax; it++){        
        // replicate flux calculation
        computeTimestepArray(solutionOld, dt);
        
        // updateMassFlows(solutionOld); for the moment no need
        
        // time-integration steps (runge-kutta in general)
        for (auto &integrationCoeff: timeIntegrationCoeffs){
            
            fluxResiduals = computeFluxResiduals(solutionOld, it);
            
            updateSolution(solutionOld, solutionNew, fluxResiduals, integrationCoeff, dt);

        }

        _conservativeVars = solutionNew;
        printInfoResiduals(fluxResiduals, it);
    }
}

void CEulerSolver::printInfoResiduals(FlowSolution &residuals, unsigned long int it) const {
    if (it == 0) {printHeader();}
    auto logRes = computeLogResidualNorm(residuals);

    std::cout << "|" << std::setw(14) << std::setfill(' ') << std::left << it+1 << "|"
                << std::setw(14) << std::left << 0.0 << "|"
                << std::setw(14) << std::right << logRes[0] << "|"
                << std::setw(14) << std::right << logRes[1] << "|"
                << std::setw(14) << std::right << logRes[2] << "|"
                << std::setw(14) << std::right << logRes[3] << "|"
                << std::setw(14) << std::right << logRes[4] << "|"
                << std::endl;
}

StateVector CEulerSolver::computeLogResidualNorm(const FlowSolution &residuals) const {
    StateVector logResidualNorm{};
    for (int i=0; i<5; i++){
        logResidualNorm[i] = std::log10(residuals.norm(i) / (_nPointsI*_nPointsJ*_nPointsK));
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

    for (int i=0; i<_nPointsI; i++) {
        for (int j=0; j<_nPointsJ; j++){
            for (int k=0; k<_nPointsK; k++){
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
    
    // for (auto &id: bcIndices){
    //     auto surfaces = _mesh.getBoundarySurface(id);
    //     // we need slices of the arrays
    //     std::cout << "Must take slices of the arrays";
    // }
    
}

FlowSolution CEulerSolver::computeFluxResiduals(const FlowSolution solution, unsigned long int it) const {
    FlowSolution residuals(_nPointsI, _nPointsJ, _nPointsK);
    computeAdvectionResiduals(FluxDirection::I, solution, it, residuals);
    computeAdvectionResiduals(FluxDirection::J, solution, it, residuals);
    computeAdvectionResiduals(FluxDirection::K, solution, it, residuals);
    return residuals;
}

void CEulerSolver::computeAdvectionResiduals(FluxDirection direction, const FlowSolution solution, unsigned long int itCounter, FlowSolution &residuals) const {
    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);
    const Matrix3D<Vector3D>& midPoints = _mesh.getMidPoints(direction);
    StateVector Uleft{}, Uright{}, Uleftleft {}, Urightright {}, flux {};
    Vector3D surface {};
    // CAdvectionScheme advectionScheme();

    auto ni = surfaces.sizeI(); 
    auto nj = surfaces.sizeJ(); 
    auto nk = surfaces.sizeK();

    for (auto iFace = 0UL; iFace < ni; ++iFace) {
        for (auto jFace = 0UL; jFace < nj; ++jFace) {
            for (auto kFace = 0UL; kFace < nk; ++kFace) {
                auto dirFace = 0UL;
                auto stopFace = 0UL;
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
                    // boundary flux calculation at the beginning
                } else if (dirFace == stopFace) {
                    // boundary flux calculation at the end
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

void CEulerSolver::updateSolution(FlowSolution &solOld, FlowSolution &solNew, const FlowSolution &residuals, const FloatType &integrationCoeff, const Matrix3D<FloatType> &dt){

    for (size_t i = 0; i < _nPointsI; i++) {
        for (size_t j = 0; j < _nPointsJ; j++) {
            for (size_t k = 0; k < _nPointsK; k++) {
                solNew.set(i, j, k, solOld.at(i, j, k) - (residuals.at(i, j, k) * integrationCoeff * dt(i, j, k) / _mesh.getVolume(i, j, k)));
            }
        }
    }
    solOld = solNew;
}

