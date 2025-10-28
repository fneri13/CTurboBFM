#include "CBoundaryConditionBase.hpp"
#include "commonFunctions.hpp"

StateVector CBoundaryConditionBase::computeSubsonicsInletFlux(StateVector internalConservative, Vector3D surface, 
    FloatType totPressure, FloatType totTemperature, Vector3D flowDirection){

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
    FloatType normalMachBound = velocityBoundMag / soundSpeedBound;
    FloatType pressureBound = _fluid.computeStaticPressure_pt_M(totPressure, normalMachBound);
    FloatType temperatureBound = _fluid.computeStaticTemperature_Tt_M(totTemperature, normalMachBound);
    FloatType densityBound = _fluid.computeDensity_p_T(pressureBound, temperatureBound);
    FloatType energyBound = _fluid.computeStaticEnergy_p_rho(pressureBound, densityBound);
    Vector3D velocityBound = flowDirection * velocityBoundMag;
    FloatType totEnergyBound = energyBound + 0.5 * velocityBound.dot(velocityBound);

    // compute boundary flux
    StateVector primitiveBoundary({densityBound, velocityBound.x(), velocityBound.y(), velocityBound.z(), totEnergyBound});
    StateVector flux = computeEulerFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;

}