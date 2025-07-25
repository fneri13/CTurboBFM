#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Thollet model source terms computation
class CSourceBFMNeri : public CSourceBFMBase {
public:

    CSourceBFMNeri(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMNeri() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

private:
    
};
