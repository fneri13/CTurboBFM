#include "CBoundaryConditionOutlet.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionOutlet::CBoundaryConditionOutlet(const Config &config, const CMesh &mesh, CFluidBase &fluid, BoundaryIndices boundIndex, std::vector<FloatType> bcValues)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = bcValues;
    }


StateVector CBoundaryConditionOutlet::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) {
    auto primitive = getEulerPrimitiveFromConservative(internalConservative);
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    auto density = primitive[0];
    auto pressure = _fluid.computePressure_rho_u_et(density, velocity, primitive[4]);
    auto soundSpeed = _fluid.computeSoundSpeed_p_rho(pressure, primitive[0]);

    if (velocity.magnitude() >= soundSpeed) {
        auto flux = computeEulerFluxFromPrimitive(primitive, surface, _fluid);
        return flux;
    }
    else {
        FloatType outletPressure = _boundaryValues[0];
        FloatType pressureBoundary = _config.computeRampedOutletPressure(iterCounter, outletPressure);
        FloatType densityBoundary = pressureBoundary * density / pressure;
        Vector3D velocityBoundary = velocity;
        FloatType energyBoundary = _fluid.computeStaticEnergy_p_rho(pressureBoundary, densityBoundary);
        FloatType totEnergyBoundary = energyBoundary + 0.5 * velocityBoundary.dot(velocityBoundary);
        StateVector primitiveBoundary({densityBoundary, velocityBoundary.x(), velocityBoundary.y(), velocityBoundary.z(), totEnergyBoundary});
        auto flux = computeEulerFluxFromPrimitive(primitiveBoundary, surface, _fluid);
        return flux;
    }
}