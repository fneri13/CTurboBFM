#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Thollet model source terms computation
class CSourceBFMLamprakis : public CSourceBFMBase {
public:

    CSourceBFMLamprakis(const Config &config, const CFluidBase &fluid, const CMesh &mesh, std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance);

    virtual ~CSourceBFMLamprakis() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    FloatType computeLoadingForceIntensity(const StateVector& primitive) const;

    /** Compute the flow state at a given point, needed for the computation of source terms
     * 
     * @param i x index
     * @param j y index
     * @param k z index
     * @param primitive primitive variables
    */
    void computeFlowState(size_t i, size_t j, size_t k, const StateVector& primitive);

    Vector3D computeLoadingVersor(const FloatType& metalAngle, const FloatType& leanAngle, const FloatType& gasPathAngle) const;

private:

    // std::map<TurboPerformance, std::vector<FloatType>> &_turboPerformance;
    FloatType _metalAngle, _leanAngle, _gasPathAngle, _d_metalAngle_dm;
    FloatType _velMeridionalMag, _betaFlow, _relativeVelocityMag;
    Vector3D _loadingVersor_cyl, _loadingVersor_cart;
    
};
