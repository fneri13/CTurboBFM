#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Gong model source terms computation
class CSourceBFMGong : public CSourceBFMBase {
public:

    CSourceBFMGong(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMGong() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    FloatType computeKnCoefficient(const FloatType normalCamberAxial, const FloatType normalCamberRadial, const FloatType normalCamberTangential) const;

private:

    // coefficients of Gong model, that have been found for the Nasa Rotor 35
    FloatType _Kn1 = 4.2;
    FloatType _Kn2 = -3.3;
    FloatType _Kp = 0.04;
    
};
