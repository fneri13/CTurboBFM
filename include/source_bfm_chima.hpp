#pragma once

#include "source_bfm_base.hpp"

// Base class handling BFM Thollet model source terms computation
class SourceBFMChima : public SourceBFMBase {
public:

    SourceBFMChima(const Config &config, const FluidBase &fluid, const Mesh &mesh, std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance);

    virtual ~SourceBFMChima() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, FlowSolution &conservativeVars) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

private:

    Vector3D _viscousForceCylindrical = {0, 0, 0};
    FloatType _tangentialForce;
    InputTable _inputTable;
    std::map<TurboPerformance, std::vector<FloatType>> &_turboPerformance;
    FloatType _scalingTurning=1.0, _scalingLoss = 1.0;
    size_t _trailingEdgeIndex = 0;
    
};
