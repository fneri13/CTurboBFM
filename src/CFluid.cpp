#include "../include/CFluid.hpp"
#include <cmath>

CFluid::CFluid(FloatType gamma, FloatType R)
    : _gamma(gamma), _R(R) {
    _cp = (_gamma * _R) / (_gamma - 1);
    _cv = _cp / _gamma;
}

FloatType CFluid::computeStaticEnergy_p_rho(FloatType p, FloatType rho) const {
    return p / ((_gamma - 1) * rho);
}

FloatType CFluid::computePressure_rho_e(FloatType rho, FloatType e) const {
    return (_gamma - 1) * rho * e;
}

FloatType CFluid::computeSoundSpeed_p_rho(FloatType p, FloatType rho) const {
    return std::sqrt(_gamma * p / rho);
}

FloatType CFluid::computeStaticEnergy_u_et(const Vector3D& vel, FloatType et) const {
    FloatType velMag2 = std::pow(vel.magnitude(), 2);
    return et - 0.5 * velMag2;
}

FloatType CFluid::computeStaticEnergy_u_et(FloatType velMag, FloatType et) const {
    return et - 0.5 * velMag * velMag;
}

FloatType CFluid::computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    FloatType p = computePressure_rho_e(rho, e);
    return computeSoundSpeed_p_rho(p, rho);
}

FloatType CFluid::computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    return computePressure_rho_e(rho, e);
}

FloatType CFluid::computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    FloatType M = computeMachNumber_rho_u_et(rho, u, et);
    return p * std::pow(1.0 + (_gamma - 1.0) / 2.0 * M * M, _gamma / (_gamma - 1));
}

FloatType CFluid::computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType T = computeTemperature_rho_u_et(rho, u, et);
    FloatType M = computeMachNumber_rho_u_et(rho, u, et);
    return T * (1.0 + (_gamma - 1.0) / 2.0 * M * M);
}

FloatType CFluid::computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    return p / (rho * _R);
}

FloatType CFluid::computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    FloatType p = computePressure_rho_e(rho, e);
    return et + p / rho;
}

FloatType CFluid::computeStaticPressure_pt_M(FloatType pt, FloatType M) const {
    return pt * std::pow(1.0 + (_gamma - 1.0) / 2.0 * M * M, -_gamma / (_gamma - 1.0));
}

FloatType CFluid::computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const {
    return Tt / (1 + (_gamma - 1) / 2 * M * M);
}

FloatType CFluid::computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    FloatType s = computeEntropy_p_rho(p, rho);
    return s;
}

FloatType CFluid::computeDensity_p_T(FloatType p, FloatType T) const {
    return p / (_R * T);
}

FloatType CFluid::computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType a = computeSoundSpeed_rho_u_et(rho, u, et);
    FloatType umag = u.magnitude();
    return umag / a;
}

FloatType CFluid::computeTotalPressure_p_M(FloatType pressure, FloatType mach) const {
    return pressure * std::pow(1.0 + (_gamma - 1.0) / 2.0 * mach * mach, _gamma / (_gamma - 1));
}

FloatType CFluid::computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const {
    return temperature * (1.0 + (_gamma - 1.0) / 2.0 * mach * mach);
}

FloatType CFluid::computeEntropy_p_rho(FloatType pressure, FloatType density) const {
    return pressure / std::pow(density, _gamma);
}

void CFluid::computeInitFields(FloatType initMach, FloatType initTemperature, FloatType initPressure, Vector3D initDirection, FloatType &density, Vector3D &velocity, FloatType &totEnergy){
    FloatType soundSpeed = std::sqrt(_gamma*_R*initTemperature);
    velocity = (initDirection / initDirection.magnitude()) * soundSpeed * initMach;
    density = computeDensity_p_T(initPressure, initTemperature);
    FloatType energy = computeStaticEnergy_p_rho(initPressure, density);
    totEnergy = energy + 0.5 * pow(velocity.magnitude(), 2);
}

FloatType CFluid::computePressure_primitive(StateVector primitive) const {
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    auto pressure = computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    return pressure;
}