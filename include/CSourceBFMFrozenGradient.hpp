#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Thollet model source terms computation
class CSourceBFMFrozenGradient : public CSourceBFMBase {
public:

    CSourceBFMFrozenGradient(const Config &config, const CFluidBase &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMFrozenGradient() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    /** Update the state variables used to compute the source terms. Reconstruct the forces starting from the gradients given in input*/
    void updateState(size_t i, size_t j, size_t k, const StateVector& primitive);

    Vector3D _viscousForceCartesian, _viscousForceCylindrical, _inviscidForceCylindrical, _inviscidForceCartesian;
    
};
