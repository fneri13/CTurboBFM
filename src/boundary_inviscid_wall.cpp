#include "boundary_inviscid_wall.hpp"
#include "math_utils.hpp"


StateVector BoundaryInviscidWall::computeBoundaryFlux(
            const StateVector& internalConservative, 
            const Vector3D& surface, 
            const Vector3D& midPoint, 
            const std::array<size_t, 3>& indices, 
            const FlowSolution& flowSolution, 
            const size_t& iterCounter) {
                
    StateVector primitive = getPrimitiveVariablesFromConservative(internalConservative);
    Vector3D velocity({primitive[1], primitive[2], primitive[3]});
    FloatType pressure = _fluid.computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    StateVector flux {{0,0,0,0,0}};
    auto Sdir = surface / surface.magnitude();
    flux[0] = 0.0;
    flux[1] = pressure * Sdir.x();
    flux[2] = pressure * Sdir.y();
    flux[3] = pressure * Sdir.z();
    flux[4] = 0.0;
    return flux;
}