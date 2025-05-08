#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Thollet model source terms computation
class CSourceBFMThollet : public CSourceBFMBase {
public:

    CSourceBFMThollet(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMThollet() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive);
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive);

    FloatType computeCompressibilityCorrection(const Vector3D& relativeVelocityCylindric, const StateVector& primitive);
    
};
