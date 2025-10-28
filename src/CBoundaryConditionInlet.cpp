#include "CBoundaryConditionInlet.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionInlet::CBoundaryConditionInlet(const Config &config, const CMesh &mesh, CFluidBase &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = inletValues;
        _referenceFrame = config.getInletReferenceFrame();
    }


StateVector CBoundaryConditionInlet::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) {
    
    // gather boundary values from config
    FloatType totPressure = _boundaryValues[0];
    FloatType totTemperature = _boundaryValues[1];
    Vector3D flowDirection = {_boundaryValues[2], _boundaryValues[3], _boundaryValues[4]};

    if (_referenceFrame == ReferenceFrame::CYLINDRICAL){
        FloatType theta = _mesh.getTheta(indices[0], indices[1], indices[2]);
        flowDirection = computeCartesianVectorFromCylindrical(flowDirection, theta);
    }

    if (flowDirection.x() == 1.0 && flowDirection.y() == 1.0 && flowDirection.z() == 1.0) {
        // in this case the direction is normal to the surface
        flowDirection = -surface;
    }

    flowDirection /= flowDirection.magnitude();

    StateVector flux = computeSubsonicsInletFlux(internalConservative, surface, totPressure, totTemperature, flowDirection);
    
    return flux;
}
