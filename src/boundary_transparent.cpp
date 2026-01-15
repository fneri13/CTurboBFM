#include "boundary_transparent.hpp"
#include "math_utils.hpp"


StateVector BoundaryTransparent::computeBoundaryFlux(
            const StateVector& internalConservative, 
            const Vector3D& surface, 
            const Vector3D& midPoint, 
            const std::array<size_t, 3>& indices, 
            const FlowSolution& flowSolution, 
            const size_t& iterCounter) {
                
    StateVector Ul = internalConservative;
    StateVector Ur = internalConservative;
    StateVector flux = _advScheme.computeFlux(Ul, Ul, Ur, Ur, surface);
    return flux;
}