#include "CBoundaryConditionWedge.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionWedge::CBoundaryConditionWedge(const Config &config, const CMesh &mesh, CFluidBase &fluid, BoundaryIndices boundIndex)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
    }


    StateVector CBoundaryConditionWedge::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) {
        StateVector primitive = getEulerPrimitiveFromConservative(internalConservative);
        
        FloatType yCoord = midPoint.y();
        FloatType zCoord = midPoint.z();
        FloatType thetaBoundary = std::atan2(zCoord, yCoord);
        
        Vector3D velocityCentralPlane = Vector3D(primitive[1], primitive[2], primitive[3]);
        Vector3D velocityBoundary = rotateVectorAlongXAxis(velocityCentralPlane, thetaBoundary);

        StateVector primitiveBoundary({primitive[0], velocityBoundary.x(), velocityBoundary.y(), velocityBoundary.z(), primitive[4]});
        StateVector flux = computeEulerFluxFromPrimitive(primitiveBoundary, surface, _fluid);
        return flux;
    }