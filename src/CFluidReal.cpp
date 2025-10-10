#include "CFluidReal.hpp"
#include <cmath>

CFluidReal::CFluidReal(std::string fluidName) : _fluidName(fluidName) {}

FloatType CFluidReal::computeStaticEnergy_p_rho(FloatType p, FloatType rho) const {
    return CoolProp::PropsSI("U", "P", p, "D", rho, _fluidName);
}

FloatType CFluidReal::computePressure_rho_e(FloatType rho, FloatType e) const {
    return CoolProp::PropsSI("P", "D", rho, "U", e, _fluidName); 
}


FloatType CFluidReal::computePressure_rho_T(FloatType rho, FloatType Temp) const {
    return CoolProp::PropsSI("P", "D", rho, "T", Temp, _fluidName); 
}

FloatType CFluidReal::computeSoundSpeed_p_rho(FloatType p, FloatType rho) const {
    return CoolProp::PropsSI("A", "P", p, "D", rho, _fluidName);
}

FloatType CFluidReal::computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u.magnitude(), et);
    return CoolProp::PropsSI("A", "D", rho, "U", e, _fluidName);
}

FloatType CFluidReal::computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u.magnitude(), et);
    return CoolProp::PropsSI("P", "D", rho, "U", e, _fluidName);
}

FloatType CFluidReal::computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    FloatType ht = et + p / rho;
    FloatType s = CoolProp::PropsSI("S", "D", rho, "P", p, _fluidName);
    return CoolProp::PropsSI("P", "H", ht, "S", s, _fluidName);
}

FloatType CFluidReal::computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    FloatType ht = et + p / rho;
    FloatType s = CoolProp::PropsSI("S", "D", rho, "P", p, _fluidName);
    return CoolProp::PropsSI("T", "H", ht, "S", s, _fluidName);
}

FloatType CFluidReal::computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    return CoolProp::PropsSI("T", "P", p, "D", rho, _fluidName);
}

FloatType CFluidReal::computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u.magnitude(), et);
    FloatType p = computePressure_rho_e(rho, e);
    return et + p / rho;
}

FloatType CFluidReal::computeStaticPressure_pt_M(FloatType pt, FloatType M) const {
    return pt/1.1; // need to be written properly
}

FloatType CFluidReal::computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const {
    return Tt/1.1; // need to be written properly
}

FloatType CFluidReal::computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    return CoolProp::PropsSI("S", "D", rho, "P", p, _fluidName);
}

FloatType CFluidReal::computeDensity_p_T(FloatType p, FloatType T) const {
    return CoolProp::PropsSI("D", "P", p, "T", T, _fluidName);
}

FloatType CFluidReal::computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType a = computeSoundSpeed_rho_u_et(rho, u, et);
    FloatType umag = u.magnitude();
    return umag / a;
}

FloatType CFluidReal::computeTotalPressure_p_M(FloatType pressure, FloatType mach) const {
    return pressure * 1.1; // need to be written properly
}

FloatType CFluidReal::computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const {
    return temperature * 1.1; // need to be written properly
}

FloatType CFluidReal::computeEntropy_p_rho(FloatType pressure, FloatType density) const {
    return CoolProp::PropsSI("S", "P", pressure, "D", density, _fluidName);
}

void CFluidReal::computeInitFields(FloatType initMach, FloatType initTemperature, FloatType initPressure, Vector3D flowDirection, FloatType &density, Vector3D &velocity, FloatType &totEnergy){
    FloatType soundSpeed = CoolProp::PropsSI("A", "P", initPressure, "T", initTemperature, _fluidName);
    velocity = (flowDirection / flowDirection.magnitude()) * soundSpeed * initMach;
    density = computeDensity_p_T(initPressure, initTemperature);
    FloatType energy = computeStaticEnergy_p_rho(initPressure, density);
    totEnergy = energy + 0.5 * pow(velocity.magnitude(), 2);
}

FloatType CFluidReal::computePressure_primitive(StateVector primitive) const {
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    auto pressure = computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    return pressure;
}

FloatType CFluidReal::computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const{
    return 1.0; // need to be written properly
}

FloatType CFluidReal::computeEntropy_p_T(FloatType pressure, FloatType temperature) const{
    return CoolProp::PropsSI("S", "P", pressure, "T", temperature, _fluidName);
}