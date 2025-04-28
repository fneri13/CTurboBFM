#pragma once
#include "types.hpp"

// Class for ideal gases
class CFluid {
public:
    CFluid(FloatType gamma, FloatType R);

    FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const;
    
    FloatType computePressure_rho_e(FloatType rho, FloatType e) const;
    
    FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const;

    FloatType computeStaticEnergy_u_et(const Vector3D& vel, FloatType et) const;
    
    FloatType computeStaticEnergy_u_et(FloatType velMag, FloatType et) const;

    FloatType computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;
    
    FloatType computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;
    
    FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;
    
    FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;

    FloatType computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;
    
    FloatType computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;

    FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const;
        
    FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const;
    
    FloatType computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;

    FloatType computeDensity_p_T(FloatType p, FloatType T) const;

    FloatType computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const;
    
    FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const;
    
    FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const;

    FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const;

    void computeInitFields(FloatType initMach, FloatType initTemperature, FloatType initPressure, Vector3D initDirection, FloatType &density, Vector3D &velocity, FloatType &totEnergy);

    FloatType computePressure_primitive(StateVector primitive) const;

    FloatType getGamma() const { return _gamma; }

private:
    FloatType _gamma;
    FloatType _R;
    FloatType _cp;
    FloatType _cv;

};
