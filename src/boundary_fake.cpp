#include "boundary_fake.hpp"
#include "math_utils.hpp"

StateVector BoundaryFake::computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) {
            
    return StateVector({0, 0, 0, 0, 0});
}