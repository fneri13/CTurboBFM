#pragma once

#include "source_bfm_base.hpp"

// Base class handling BFM Thollet model source terms computation
class SourceBFMThollet : public SourceBFMBase {
public:

    SourceBFMThollet(const Config &config, const FluidBase &fluid, const Mesh &mesh) 
        : SourceBFMBase(config, fluid, mesh) {}

    virtual ~SourceBFMThollet() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, 
                Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, FlowSolution &conservativeVars);

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    FloatType computeCompressibilityCorrection(const Vector3D& relativeVelocityCylindric, const StateVector& primitive);

private:
    FloatType _Kn = _config.getHallTholletCoefficient_KN();
    FloatType _Kf = _config.getHallTholletCoefficient_KF();
    FloatType _Kd = _config.getHallTholletCoefficient_KD();
    std::string _offDesignActive = _config.getHallTholletOffDesignActive();
    
};
