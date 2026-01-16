#include "boundary_base.hpp"
#include "math_utils.hpp"

StateVector BoundaryBase::computeSubsonicInletFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const FloatType& totPressureBoundary, 
    const FloatType& totTemperatureBoundary, 
    const Vector3D& flowDirection){

    // properties of internal point
    StateVector primitive = getPrimitiveVariablesFromConservative(internalConservative);
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
    FloatType pressureBound = _fluid.computeStaticPressure_pt_M(totPressureBoundary, normalMachBound);
    FloatType temperatureBound = _fluid.computeStaticTemperature_Tt_M(totTemperatureBoundary, normalMachBound);
    FloatType densityBound = _fluid.computeDensity_p_T(pressureBound, temperatureBound);
    FloatType energyBound = _fluid.computeStaticEnergy_p_rho(pressureBound, densityBound);
    Vector3D velocityBound = flowDirection * velocityBoundMag;
    FloatType totEnergyBound = energyBound + 0.5 * velocityBound.dot(velocityBound);

    // compute boundary flux
    StateVector primitiveBoundary({
        densityBound, 
        velocityBound.x(), 
        velocityBound.y(), 
        velocityBound.z(), 
        totEnergyBound});
        
    StateVector flux = computeAdvectionFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;

}

StateVector BoundaryBase::computeOutletFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const FloatType& boundaryPressure,
    const size_t& iterCounter){

    auto primitive = getPrimitiveVariablesFromConservative(internalConservative);
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    auto density = primitive[0];
    auto pressure = _fluid.computePressure_rho_u_et(density, velocity, primitive[4]);
    auto soundSpeed = _fluid.computeSoundSpeed_p_rho(pressure, primitive[0]);
    
    if (velocity.magnitude() >= soundSpeed) {
        auto flux = computeAdvectionFluxFromPrimitive(primitive, surface, _fluid);
        return flux;
    }
    else {
        FloatType pressureBoundary = _config.computeRampedOutletPressure(iterCounter, boundaryPressure);
        FloatType densityBoundary = pressureBoundary * density / pressure;
        Vector3D velocityBoundary = velocity;
        FloatType energyBoundary = _fluid.computeStaticEnergy_p_rho(pressureBoundary, densityBoundary);
        FloatType totEnergyBoundary = energyBoundary + 0.5 * velocityBoundary.dot(velocityBoundary);
        StateVector primitiveBoundary({
            densityBoundary, 
            velocityBoundary.x(), 
            velocityBoundary.y(), 
            velocityBoundary.z(), 
            totEnergyBoundary});
        auto flux = computeAdvectionFluxFromPrimitive(primitiveBoundary, surface, _fluid);
        return flux;
    }

}

