#include "CBoundaryConditionPeriodic.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionPeriodic::CBoundaryConditionPeriodic(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = inletValues;
    }


StateVector CBoundaryConditionPeriodic::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution) {
    // properties of internal point
    StateVector primitive = getEulerPrimitiveFromConservative(internalConservative);
    Vector3D velocityInt({primitive[1], primitive[2], primitive[3]});
    FloatType soundSpeedInt = _fluid.computeSoundSpeed_rho_u_et(primitive[0], velocityInt, primitive[4]);
    FloatType totEnthalpyInt = _fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocityInt, primitive[4]);
    FloatType Jm = - velocityInt.magnitude() + 2*soundSpeedInt / (_fluid.getGamma() - 1);

    // Solve the quadratic equation for the speed of sound
    FloatType alpha = 1.0 / (_fluid.getGamma() - 1.0) + 2.0 / std::pow((_fluid.getGamma() - 1.0), 2);
    FloatType beta = -2.0 * Jm / (_fluid.getGamma() - 1.0);
    FloatType zeta = 0.5 * Jm * Jm - totEnthalpyInt;
    FloatType soundSpeedBound = std::max((-beta + std::sqrt(beta*beta - 4.0*alpha*zeta))/2.0/alpha,
                                         (-beta - std::sqrt(beta*beta - 4.0*alpha*zeta))/2.0/alpha);

    // reconstruct the boundary state                                     
    FloatType velocityBoundMag = 2.0*soundSpeedBound / (_fluid.getGamma() - 1.0) - Jm;
    Vector3D flowDirection = {_boundaryValues[2], _boundaryValues[3], _boundaryValues[4]};
    flowDirection /= flowDirection.magnitude();
    FloatType normalMachBound = velocityBoundMag / soundSpeedBound;
    FloatType pressureBound = _fluid.computeStaticPressure_pt_M(_boundaryValues[0], normalMachBound);
    FloatType temperatureBound = _fluid.computeStaticTemperature_Tt_M(_boundaryValues[1], normalMachBound);
    FloatType densityBound = _fluid.computeDensity_p_T(pressureBound, temperatureBound);
    FloatType energyBoubd = _fluid.computeStaticEnergy_p_rho(pressureBound, densityBound);
    Vector3D velocityBound = flowDirection * velocityBoundMag;
    FloatType totEnergyBound = energyBoubd + 0.5 * velocityBound.dot(velocityBound);

    // compute boundary flux
    StateVector primitiveBoundary({densityBound, velocityBound.x(), velocityBound.y(), velocityBound.z(), totEnergyBound});
    auto flux = computeEulerFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;
}