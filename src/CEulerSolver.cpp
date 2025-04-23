#include "../include/CEulerSolver.hpp"
#include "../include/commonFunctions.hpp"
#include "../include/types.hpp"
#include <iostream> // Optional, for logging/debugging
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
                _conservativeVars.rho(i,j,k) = density;
                _conservativeVars.rhoU(i,j,k) = density * velocity.x();
                _conservativeVars.rhoV(i,j,k) = density * velocity.y();
                _conservativeVars.rhoW(i,j,k) = density * velocity.z();
                _conservativeVars.rhoE(i,j,k) = density * totEnergy;
            }
        }
    }

}

void CEulerSolver::solve(){
    unsigned long int nIterMax = _config.getMaxIterations();
    FloatType physicalTime {0.0};
    Matrix3D<FloatType> dt(_nPointsI, _nPointsJ, _nPointsK);

    for (unsigned long int it=0; it<nIterMax; it++){
        std::cout << "Iteration number " << it+1 << std::endl;
        EulerConservativeVariables solutionOld = _conservativeVars;
        computeTimestepArray(solutionOld, dt);
        updateMassFlows(solutionOld);


    }
}

Matrix3D<FloatType> CEulerSolver::computeTimestepArray(const EulerConservativeVariables &solution, Matrix3D<FloatType> &timestep){
    FloatType cflMax = _config.getCFL();
    Vector3D iEdge, jEdge, kEdge;
    Vector3D iDir, jDir, kDir;
    Vector3D velocity;
    std::array<FloatType, 5> primitive{0.0}, conservative{0.0};
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


void CEulerSolver::updateMassFlows(const EulerConservativeVariables &solution){
    std::array<BoundaryIndices, 6> bcIndices {BoundaryIndices::I_START,
                                              BoundaryIndices::I_END,
                                              BoundaryIndices::J_START,
                                              BoundaryIndices::J_END,
                                              BoundaryIndices::K_START,
                                              BoundaryIndices::K_END};
    
    for (auto &id: bcIndices){
        auto surfaces = _mesh.getBoundarySurface(id);
        // we need slices of the arrays
        std::cout << "Must take slices of the arrays";
    }
    
}