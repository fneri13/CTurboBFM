#include "boundary_outlet_radial_equilibrium.hpp"
#include "math_utils.hpp"


StateVector BoundaryOutletRadialEquilibrium::computeBoundaryFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const Vector3D& midPoint, 
    const std::array<size_t, 3>& indices, 
    const FlowSolution& flowSolution, 
    const size_t& iterCounter) {

    FloatType outletPressure = _radialPressureProfile[indices[1]];
    auto flux = computeOutletFlux(internalConservative, surface, outletPressure, iterCounter);
    return flux;
}
