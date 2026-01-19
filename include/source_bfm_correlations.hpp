#pragma once

#include "source_bfm_base.hpp"

class SourceBFMCorrelations : public SourceBFMBase {
public:

    SourceBFMCorrelations(const Config &config, const FluidBase &fluid, const Mesh &mesh);

    virtual ~SourceBFMCorrelations() = default;  

protected:

    virtual StateVector computeBodyForceSource(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce, 
        Matrix3D<Vector3D> &viscousForce, 
        FlowSolution &conservativeVars);

    StateVector computeInviscidComponent(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce);
    
    /** @brief not implemented yet --> zero viscous force for now*/
    StateVector computeViscousComponent(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &viscousForce);

private:

    void computeCorrelationParameters(
        size_t i, 
        size_t j, 
        size_t k, 
        FlowSolution &conservativeSolution);

private:
    size_t _leadingEdgeIdx;
    size_t _trailingEdgeIdx;
    FloatType _solidity;
    Vector3D _relVelInlet, _relVelOutlet;
    FloatType _diffusionFactor; // Lieblein definition
    FloatType _inletFlowAngleDeg;
    FloatType _bladeMeridionalLength;
    FloatType _camberAngleDegAbs;
    FloatType _incidenceZeroTenThkStar;
    
};
