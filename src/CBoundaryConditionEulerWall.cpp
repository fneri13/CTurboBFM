#include "../include/CBoundaryConditionEulerWall.hpp"
#include "../include/commonFunctions.hpp"

CBoundaryConditionEulerWall::CBoundaryConditionEulerWall(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        std::cout << "Euler wall boundary condition initialized\n";
    }


    StateVector CBoundaryConditionEulerWall::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint) {
        auto primitive = getEulerPrimitiveFromConservative(internalConservative);
        auto pressure = _fluid.computePressure_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);
        StateVector flux{};
        auto Sdir = surface / surface.magnitude();
        flux[0] = 0.0;
        flux[1] = pressure * Sdir.x();
        flux[2] = pressure * Sdir.y();
        flux[3] = pressure * Sdir.z();
        flux[4] = 0.0;
        return flux;
    }