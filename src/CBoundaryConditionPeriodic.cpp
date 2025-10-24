#include "CBoundaryConditionPeriodic.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionPeriodic::CBoundaryConditionPeriodic(const Config &config, const CMesh &mesh, CFluidBase &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = inletValues;
        _periodicityAngle = _boundaryValues.at(0); // just an alias for it
    }


StateVector CBoundaryConditionPeriodic::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &solution, const size_t iterCounter) {
    // the real way to do it would be to use the numerical flux schemes using internal points, and then using the logic described in blazek. This case
    // introduces a bit of differences because the flux is computed directly from the state on the boundary, so the fluxes dont compensate exactly
    return StateVector({0, 0, 0, 0, 0});

    // // properties of the periodic node
    // size_t periodicIdx = 0;
    // FloatType thetaRotation = 0;
    // if (indices[2] == 0) {
    //     periodicIdx = _mesh.getNumberPointsK() - 1;
    //     thetaRotation = - _periodicityAngle;
    // } 
    // else {
    //     periodicIdx = 0;
    //     thetaRotation = _periodicityAngle;
    // }

    // StateVector conservativePeriodic = solution.at(indices[0], indices[1], periodicIdx);
    // StateVector primitivePeriodic = getEulerPrimitiveFromConservative(conservativePeriodic);
    // Vector3D velocityPeriodic({primitivePeriodic[1], primitivePeriodic[2], primitivePeriodic[3]});

    
    // // rotate that vector to be on the boundary where the flux is computed
    // Vector3D velocityBoundaryNew = rotateVectorAlongXAxis(velocityPeriodic, thetaRotation);

    // StateVector primitiveBoundaryNew({primitivePeriodic[0], velocityBoundaryNew.x(), velocityBoundaryNew.y(), velocityBoundaryNew.z(), primitivePeriodic[4]});
    // StateVector flux = computeEulerFluxFromPrimitive(primitiveBoundaryNew, surface, _fluid);
    
    // return flux;
}