#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Correlations model source terms computation
class CSourceBFMCorrelations : public CSourceBFMBase {
public:

    CSourceBFMCorrelations(const Config &config, const CFluidBase &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMCorrelations() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, 
                Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, FlowSolution &conservativeVars);

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

private:
    FloatType _incidenceZeroTenThkStar;
    
};
