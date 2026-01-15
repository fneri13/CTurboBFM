#include "boundary_outlet_supersonic.hpp"
#include "math_utils.hpp"

StateVector BoundaryOutletSupersonic::computeBoundaryFlux(
            const StateVector& internalConservative, 
            const Vector3D& surface, 
            const Vector3D& midPoint, 
            const std::array<size_t, 3>& indices, 
            const FlowSolution& flowSolution, 
            const size_t& iterCounter) {
                
    auto primitive = getEulerPrimitiveFromConservative(internalConservative);
    auto flux = computeEulerFluxFromPrimitive(primitive, surface, _fluid);
    return flux;
}