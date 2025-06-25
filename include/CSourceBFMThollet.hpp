#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Thollet model source terms computation
class CSourceBFMThollet : public CSourceBFMBase {
public:

    CSourceBFMThollet(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMThollet() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    FloatType computeCompressibilityCorrection(const Vector3D& relativeVelocityCylindric, const StateVector& primitive);

private:
    FloatType _Kn = _config.getHallTholletCoefficient_KN();
    FloatType _Kf = _config.getHallTholletCoefficient_KF();
    FloatType _Kd = _config.getHallTholletCoefficient_KD();
    
};
