#pragma once

#include "source_bfm_base.hpp"

class SourceBFMChima : public SourceBFMBase {
public:

    SourceBFMChima(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh, 
        std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance);

    virtual ~SourceBFMChima() = default;  

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
        Matrix3D<Vector3D> &inviscidForce,
        FloatType &dt,
        FlowSolution &conservativeVars);
    
    StateVector computeViscousComponent(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &viscousForce);
    
    FloatType computeDerivativeFunction(FloatType finalValue, FloatType xLocal, FloatType xFinal);

private:
    InputTable _inputTable;
    std::map<TurboPerformance, std::vector<FloatType>> &_turboPerformance;
    FloatType _scalingTurning=1.0, _scalingLoss = 1.0;
    size_t _leadingEdgeIndex = 0;
    size_t _trailingEdgeIndex = 0;
    Vector3D _viscousForceCylindrical = {0, 0, 0};
    FloatType _tangentialForce;
    
};
