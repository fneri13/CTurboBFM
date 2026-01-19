#pragma once

#include "source_bfm_base.hpp"

class SourceBFMHall : public SourceBFMBase {
public:

    SourceBFMHall(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh) 
        : SourceBFMBase(config, fluid, mesh) {}

    virtual ~SourceBFMHall() = default;  

protected:

    virtual StateVector computeBodyForceSource(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce, 
        Matrix3D<Vector3D> &viscousForce, 
        FlowSolution &conservativeVars) override;
    
};
