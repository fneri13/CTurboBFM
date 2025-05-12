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
    
    _output = std::make_unique<COutputCSV>(_config, _mesh, _conservativeVars, *_fluid, _inviscidForce, _viscousForce);

    BFM_Model bfmModel = _config.getBFMModel();
    if (bfmModel == BFM_Model::HALL) {
        _bfmSource = std::make_unique<CSourceBFMHall>(_config, *_fluid, _mesh);
    }
    else if (bfmModel == BFM_Model::HALL_THOLLET) {
        _bfmSource = std::make_unique<CSourceBFMThollet>(_config, *_fluid, _mesh);
    }
    else {
        _bfmSource = std::make_unique<CSourceBFMBase>(_config, *_fluid, _mesh);
    }
    
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

    _radialProfilePressure.resize(_nPointsJ);
    _radialProfileCoords.resize(_nPointsJ);

    _inviscidForce.resize(_nPointsI, _nPointsJ, _nPointsK);
    _viscousForce.resize(_nPointsI, _nPointsJ, _nPointsK);

    size_t nj = _mesh.getNumberPointsJ();
    for (size_t j = 0; j < nj; j++) {
        _radialProfileCoords[j] = std::sqrt(_mesh.getVertex(0,j,0).y()*_mesh.getVertex(0,j,0).y() + _mesh.getVertex(0,j,0).z()*_mesh.getVertex(0,j,0).z());
    }

}

void CEulerSolver::solve(){
    size_t nIterMax = _config.getMaxIterations();
    Matrix3D<FloatType> timestep(_nPointsI, _nPointsJ, _nPointsK);                          // place holder for time step array
    std::vector<FloatType> timeIntegrationCoeffs = _config.getTimeIntegrationCoeffs();      // time integration coefficients (runge kutta)
    FlowSolution fluxResiduals(_nPointsI, _nPointsJ, _nPointsK);                            // place holder for flux residuals
    size_t updateMassFlowsFreq = _config.getMassFlowUpdateFrequency();                      // frequency to update the mass flows at the boundaries
    size_t solutionOutputFreq = _config.getSolutionOutputFrequency();                       // frequency to output the solution
    bool turboOutput = _config.saveTurboOutput();                                           // flag to save the solution in turbo format

    // time integration
    for (size_t it=1; it<nIterMax; it++){        
        FlowSolution solutionOld = _conservativeVars;                                       // place holder for the solution at the previous timestep
        if (it%updateMassFlowsFreq == 0) {
            updateMassFlows(solutionOld);
            if (turboOutput) updateTurboPerformance(solutionOld);
        }                      
        computeTimestepArray(solutionOld, timestep);                                        // compute the physical time step
        FlowSolution tmpSol = solutionOld;                                                  // place holder for the solution at the runge-kutta step
        
        // runge-kutta steps
        for (const auto &integrationCoeff: timeIntegrationCoeffs){
            updateRadialProfiles(tmpSol);
            fluxResiduals = computeFluxResiduals(tmpSol, it);
            updateSolution(solutionOld, tmpSol, fluxResiduals, integrationCoeff, timestep);
        }
        
        // update the solution and print information
        _conservativeVars = tmpSol;
        printInfoResiduals(fluxResiduals, it);
        if (it%updateMassFlowsFreq == 0) {
            printInfoMassFlows(it);
            if (turboOutput) printInfoTurboPerformance(it);
        }
        if (it%solutionOutputFreq == 0) _output->writeSolution(); 
        if (it%solutionOutputFreq == 0) {
            writeLogResidualsToCSV();
        } 
    }
}

void CEulerSolver::printInfoResiduals(FlowSolution &residuals, size_t it) {
    if (it == 1) {printHeader();}
    auto logRes = computeLogResidualNorm(residuals);
    printLogResiduals(logRes, it);
    _logResiduals.push_back(logRes);
}


void CEulerSolver::printInfoMassFlows(size_t it) const {
    std::cout << "\nMASS FLOWS CHECK [kg/s]:\n";
    std::cout << "I_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::I_START) << std::endl;
    std::cout << "I_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::I_END) << std::endl;
    std::cout << "J_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::J_START) << std::endl;
    std::cout << "J_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::J_END) << std::endl;
    std::cout << "K_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::K_START) << std::endl;
    std::cout << "K_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndices::K_END) << std::endl << std::endl;
}

void CEulerSolver::printInfoTurboPerformance(size_t it) const {
    std::cout << "\nTURBOMACHINERY PERFORMANCE:\n";
    std::cout << "Mass Flow [kg/s]: " << std::setprecision(6) << _turboPerformance.at(TurboPerformance::MASS_FLOW) << std::endl;
    std::cout << "Total Pressure Ratio [-]: " << std::setprecision(6) << _turboPerformance.at(TurboPerformance::TOTAL_PRESSURE_RATIO) << std::endl;
    std::cout << "Total Temperature Ratio [-]: " << std::setprecision(6) << _turboPerformance.at(TurboPerformance::TOTAL_TEMPERATURE_RATIO) << std::endl;
    std::cout << "Total Efficiency [-]: " << std::setprecision(6) << _turboPerformance.at(TurboPerformance::TOTAL_EFFICIENCY) << std::endl << std::endl;
}

void CEulerSolver::printLogResiduals(const StateVector &logRes, unsigned long int it) const {
    int col_width = 14;
    std::cout << std::fixed << std::setprecision(6); // Set fixed format with 6 decimals
    std::cout << "|" << std::setw(col_width) << std::setfill(' ') << std::left << it << "|"
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
    constexpr double minDouble = std::numeric_limits<FloatType>::min();

    for (int i = 0; i < 5; i++) {
        auto residualNorm = residuals.norm(i);
        if (residualNorm >= minDouble) {
            logResidualNorm[i] = std::log10(residualNorm / (_nPointsI * _nPointsJ * _nPointsK));
        } else {
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

void CEulerSolver::computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep){
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

    // if time step method is global, the transient is coherent, and dt must be the same for every cell
    if (_config.getTimeStepMethod() == TimeStepMethod::GLOBAL){
        dtMin = timestep.min();
        for (size_t i=0; i<_nPointsI; i++) {
            for (size_t j=0; j<_nPointsJ; j++){
                for (size_t k=0; k<_nPointsK; k++){
                    timestep(i,j,k) = dtMin;
                }
            }
        }
    }
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


void CEulerSolver::updateTurboPerformance(const FlowSolution&solution){
    _turboPerformance[TurboPerformance::MASS_FLOW] = 0.5 * (_massFlows[BoundaryIndices::I_START] + _massFlows[BoundaryIndices::I_END]);
    if (_config.getTopology() == Topology::AXISYMMETRIC){
        _turboPerformance[TurboPerformance::MASS_FLOW] *= 2.0 * M_PI / _mesh.getWedgeAngle();
    }
    

    std::array<BoundaryIndices, 2> bcIndices {BoundaryIndices::I_START, BoundaryIndices::I_END};
    std::vector<FloatType> totalPressure;
    std::vector<FloatType> totalTemperature;

    for (auto& bcIndex: bcIndices){
        Matrix2D<Vector3D> surface = _mesh.getMeshBoundary(bcIndex);

        size_t nj = surface.sizeI();
        size_t nk = surface.sizeJ();
        Matrix2D<FloatType> rhoUX_pt(nj, nk);
        Matrix2D<FloatType> rhoUY_pt(nj, nk);
        Matrix2D<FloatType> rhoUZ_pt(nj, nk);
        Matrix2D<FloatType> rhoUX_Tt(nj, nk);
        Matrix2D<FloatType> rhoUY_Tt(nj, nk);
        Matrix2D<FloatType> rhoUZ_Tt(nj, nk);

        for (size_t j=0; j<nj; j++){
            for (size_t k=0; k<nk; k++){
                StateVector primitive;
                if (bcIndex == BoundaryIndices::I_START){
                    primitive = getEulerPrimitiveFromConservative(solution.at(0,j,k));
                }
                else{
                    primitive = getEulerPrimitiveFromConservative(solution.at(_nPointsI-1,j,k));
                }
                FloatType rho = primitive[0];
                FloatType ux = primitive[1];
                FloatType uy = primitive[2];
                FloatType uz = primitive[3];
                FloatType et = primitive[4];
                FloatType totalPressure = _fluid->computeTotalPressure_rho_u_et(rho, {ux,uy,uz}, et);
                FloatType totalTemperature = _fluid->computeTotalTemperature_rho_u_et(rho, {ux,uy,uz}, et);

                rhoUX_pt(j,k) = rho * ux * totalPressure;
                rhoUY_pt(j,k) = rho * uy * totalPressure;
                rhoUZ_pt(j,k) = rho * uz * totalPressure;
                rhoUX_Tt(j,k) = rho * ux * totalTemperature;
                rhoUY_Tt(j,k) = rho * uy * totalTemperature;
                rhoUZ_Tt(j,k) = rho * uz * totalTemperature;
            }
        }

        totalPressure.push_back(computeSurfaceIntegral(surface, rhoUX_pt, rhoUY_pt, rhoUZ_pt) / _massFlows[bcIndex]);
        totalTemperature.push_back(computeSurfaceIntegral(surface, rhoUX_Tt, rhoUY_Tt, rhoUZ_Tt) / _massFlows[bcIndex]);
    }

    _turboPerformance[TurboPerformance::TOTAL_PRESSURE_RATIO] = totalPressure.at(1) / totalPressure.at(0);
    _turboPerformance[TurboPerformance::TOTAL_TEMPERATURE_RATIO] = totalTemperature.at(1) / totalTemperature.at(0);
    _turboPerformance[TurboPerformance::TOTAL_EFFICIENCY] = _fluid->computeTotalEfficiency_PRtt_TRt(_turboPerformance[TurboPerformance::TOTAL_PRESSURE_RATIO], _turboPerformance[TurboPerformance::TOTAL_TEMPERATURE_RATIO]);
    
}

FlowSolution CEulerSolver::computeFluxResiduals(const FlowSolution& solution, size_t it) {
    FlowSolution residuals(_nPointsI, _nPointsJ, _nPointsK); // residuals place-holder, passed by reference to below functions

    // compute residuals contribution from advection
    computeAdvectionResiduals(FluxDirection::I, solution, it, residuals);
    computeAdvectionResiduals(FluxDirection::J, solution, it, residuals);
    if (_topology==Topology::THREE_DIMENSIONAL || _topology==Topology::AXISYMMETRIC){
        computeAdvectionResiduals(FluxDirection::K, solution, it, residuals);
    }

    // compute residuals contribution from sources
    computeSourceResiduals(solution, it, residuals, _inviscidForce, _viscousForce);

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
                    flux = _boundaryConditions.at(boundaryStart)->computeBoundaryFlux(Uinternal, surface, midPoint, {iFace, jFace, kFace});
                    residuals.add(iFace, jFace, kFace, flux * surface.magnitude());
                } else if (dirFace == stopFace) {
                    Uinternal = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
                    surface = surfaces(iFace, jFace, kFace);
                    midPoint = midPoints(iFace, jFace, kFace);
                    flux = _boundaryConditions.at(boundaryEnd)->computeBoundaryFlux(Uinternal, surface, midPoint, {iFace, jFace, kFace});
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


void CEulerSolver::writeLogResidualsToCSV() const {

    std::string filename = "residuals.csv";
    std::ofstream file(filename); // open in truncate (default) mode

    if (!file.is_open()) {
        std::cerr << "Error: Could not open log residuals file: " << filename << std::endl;
        return;
    }

    // Write header
    file << "Rho,RhoU,RhoV,RhoW,RhoE\n";

    size_t size = _logResiduals.size();
    std::array<FloatType, 5> logRes;
    for (size_t i = 0; i < size; i++) {
        logRes[0] = _logResiduals[i][0];
        logRes[1] = _logResiduals[i][1];
        logRes[2] = _logResiduals[i][2];
        logRes[3] = _logResiduals[i][3];
        logRes[4] = _logResiduals[i][4];
        // Write data row
        file << logRes[0] << "," << logRes[1] << "," << logRes[2] << "," << logRes[3] << "," << logRes[4] << "\n";
    }

    file.close();
}


void CEulerSolver::updateRadialProfiles(FlowSolution &solution){
    size_t nj = _mesh.getNumberPointsJ();
    StateVector conservative, primitive;
    Vector3D velocityCart, velocityCyl;
    std::vector<FloatType> densityProfile(nj);
    std::vector<FloatType> velTangProfile(nj);
    FloatType theta;

    size_t ni = _mesh.getNumberPointsI();
    for (size_t j = 0; j < nj; j++) {
        conservative = solution.at(ni-1, j, 0);
        primitive = getEulerPrimitiveFromConservative(conservative);    
        velocityCart(0) = primitive[1];
        velocityCart(1) = primitive[2];
        velocityCart(2) = primitive[3];

        densityProfile[j] = primitive[0];
        theta = std::atan2(_mesh.getVertex(0,j,0).z(), _mesh.getVertex(0,j,0).y());
        velocityCyl = computeCylindricalVectorFromCartesian(velocityCart, theta);
        velTangProfile[j] = velocityCyl.z();
    }

    integrateRadialEquilibrium(densityProfile, velTangProfile, _radialProfileCoords, _hubStaticPressure, _radialProfilePressure);

}


void CEulerSolver::computeSourceResiduals(const FlowSolution& solution, size_t itCounter, FlowSolution &residuals, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    if (!_config.isBFMActive()){
        return ;
    }

    for (size_t i = 0; i < _nPointsI; i++) {
        for (size_t j = 0; j < _nPointsJ; j++) {
            for (size_t k = 0; k < _nPointsK; k++) {
                Vector3D blockageGradient = _mesh.getInputFieldsGradient(FieldNames::BLOCKAGE, i, j, k); 
                if (blockageGradient.magnitude() > 1E-10){
                    StateVector conservative = solution.at(i, j, k);
                    StateVector primitive = getEulerPrimitiveFromConservative(conservative);
                    StateVector source = _bfmSource->computeSource(i, j, k, primitive, inviscidForce, viscousForce);
                    residuals.subtract(i, j, k, source);
                }
            }
        }
    }
}
