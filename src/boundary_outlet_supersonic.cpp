#include "boundary_outlet_supersonic.hpp"
#include "math_utils.hpp"

StateVector BoundaryOutletSupersonic::computeBoundaryFlux(
            const StateVector& internalConservative, 
            const Vector3D& surface, 
            const Vector3D& midPoint, 
            const std::array<size_t, 3>& indices, 
            const FlowSolution& flowSolution, 
            const size_t& iterCounter) {
                
    auto primitive = getPrimitiveVariablesFromConservative(internalConservative);
    auto flux = computeAdvectionFluxFromPrimitive(primitive, surface, _fluid);
    return flux;
}