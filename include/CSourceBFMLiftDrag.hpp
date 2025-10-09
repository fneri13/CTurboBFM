#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Hall model source terms computation
class CSourceBFMLiftDrag : public CSourceBFMBase {
public:

    CSourceBFMLiftDrag(const Config &config, const CFluidBase &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMLiftDrag() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;
    
    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);

    StateVector computeInviscidComponentGongStyle(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    void computeRelativeFlowAngle(size_t i, size_t j, size_t k, const StateVector& primitive);

    FloatType _relativeFlowAngle;
};
