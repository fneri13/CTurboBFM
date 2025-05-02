#include "CBoundaryConditionEulerWall.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionEulerWall::CBoundaryConditionEulerWall(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {}


    StateVector CBoundaryConditionEulerWall::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices) {
        StateVector primitive = getEulerPrimitiveFromConservative(internalConservative);
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