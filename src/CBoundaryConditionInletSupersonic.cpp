#include "CBoundaryConditionInletSupersonic.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionInletSupersonic::CBoundaryConditionInletSupersonic(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = inletValues;
    }


StateVector CBoundaryConditionInletSupersonic::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint) {
    FloatType pressureBoundary = _boundaryValues[0];
    FloatType temperatureBoundary = _boundaryValues[1];
    Vector3D velocityBoundary({_boundaryValues[2], _boundaryValues[3], _boundaryValues[4]});
    FloatType densityBoundary = _fluid.computeDensity_p_T(pressureBoundary, temperatureBoundary);
    FloatType energyBoundary = _fluid.computeStaticEnergy_p_rho(pressureBoundary, densityBoundary);
    FloatType totEnergyBoundary = energyBoundary + 0.5 * velocityBoundary.dot(velocityBoundary);
    StateVector primitiveBoundary({densityBoundary, velocityBoundary.x(), velocityBoundary.y(), velocityBoundary.z(), totEnergyBoundary});
    auto flux = computeEulerFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;
}
