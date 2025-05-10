#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Hall model source terms computation
class CSourceBFMHall : public CSourceBFMBase {
public:

    CSourceBFMHall(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMHall() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;
    
};
