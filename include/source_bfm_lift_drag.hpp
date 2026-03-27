#pragma once

#include "source_bfm_base.hpp"

class SourceBFMLiftDrag : public SourceBFMBase {
public:

    SourceBFMLiftDrag(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh) 
        : SourceBFMBase(config, fluid, mesh) {}

    virtual ~SourceBFMLiftDrag() = default;  

protected:

    virtual StateVector computeBodyForceSource(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce, 
        Matrix3D<Vector3D> &viscousForce, 
        FlowSolution &conservativeVars,
        FloatType &dt,
        FloatType &timePhysical) override;

    StateVector computeInviscidComponent(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &viscousForce);

private:
    FloatType _solidity;

};

