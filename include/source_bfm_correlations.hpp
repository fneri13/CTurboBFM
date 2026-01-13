#pragma once

#include "source_bfm_base.hpp"

// Base class handling BFM Correlations model source terms computation
class SourceBFMCorrelations : public SourceBFMBase {
public:

    SourceBFMCorrelations(const Config &config, const FluidBase &fluid, const Mesh &mesh) : SourceBFMBase(config, fluid, mesh) {
        _leadingEdgeIdx = _config.getLeadingEdgeIndex();
        _trailingEdgeIdx = _config.getTrailingEdgeIndex();
    }

    virtual ~SourceBFMCorrelations() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, 
                Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, FlowSolution &conservativeVars);

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

private:
    FloatType _incidenceZeroTenThkStar;
    
};
