#include "boundary_outlet.hpp"
#include "math_utils.hpp"

StateVector BoundaryOutlet::computeBoundaryFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const Vector3D& midPoint, 
    const std::array<size_t, 3>& indices, 
    const FlowSolution& flowSolution, 
    const size_t& iterCounter) {
    
    FloatType outletPressure = _boundaryValues[0];
    auto flux = computeOutletFlux(internalConservative, surface, outletPressure, iterCounter);
    return flux;
}