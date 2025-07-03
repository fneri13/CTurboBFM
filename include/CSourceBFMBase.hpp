#pragma once

#include "types.hpp"
#include "CMesh.hpp"
#include "Config.hpp"
#include "CFluid.hpp"
#include "commonFunctions.hpp"

// Base class handling BFM source terms computation
class CSourceBFMBase {
public:

    CSourceBFMBase(const Config &config, const CFluid &fluid, const CMesh &mesh);

    virtual ~CSourceBFMBase() = default;  

    virtual StateVector computeSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, Matrix3D<FloatType> &deviationAngle, FloatType &timePhysical);

protected:
    const Config& _config;
    
    const CFluid& _fluid;
    
    const CMesh& _mesh;

    virtual StateVector computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive);

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce);

    Vector3D _point;
    FloatType _radius;
    FloatType _theta;
    FloatType _blockage, _bladeIsPresent;
    FloatType _omega;
    FloatType _numberBlades;
    FloatType _pitch;
    Vector3D _velocityCartesian, _velocityCylindrical;
    Vector3D _relativeVelocityCartesian, _dragVelocityCylindrical, _relativeVelocityCylindric;
    FloatType _normalCamberAxial, _normalCamberRadial, _normalCamberTangential;
    Vector3D _normalCamberCylindric;
    FloatType _deviationAngle;
    Vector3D _inviscidForceDirectionCylindrical, _inviscidForceDirectionCartesian;
    Vector3D _viscousForceDirectionCylindrical, _viscousForceDirectionCartesian;
    bool _perturbationBodyForceActive = false;

    // body force perturbation data
    FloatType _perturbationRadius = 0.0;
    Vector3D _perturbationCenter = {0.0, 0.0, 0.0};
    FloatType _perturbationScalingFactor = 1.0;
    FloatType _perturbationTimeStart = 0.0;
    FloatType _perturbationTimeDuration = 0.0;

    /** Compute the flow state at a given point, needed for the computation of source terms
     * 
     * @param i x index
     * @param j y index
     * @param k z index
     * @param primitive primitive variables
    */
    void computeFlowState(size_t i, size_t j, size_t k, const StateVector& primitive);

    /** Compute the deviation angle, using the relative velocity and the normal camber 
     * 
     * @param relativeVelocity relative velocity
     * @param normalCamber normal camber
    */
    FloatType computeDeviationAngle(Vector3D relativeVelocity, Vector3D normalCamber);

    /** Compute the inviscid force direction, using the relative velocity and the normal camber. The direction must be orthogonal to the relative velocity, 
     * and the radial component of it should reflect the radial component of the normal camber. Non linear set of equations needed here, problems possile.
     * 
     * @param relativeVelocity relative velocity
     * @param normalCamber normal camber
     */
    Vector3D computeInviscidForceDirection(const Vector3D& relativeVelocity, const Vector3D& normalCamber);

    /** Compute the tangential component of the inviscid force direction. Function used from the other method */
    FloatType computeTangentialComponent(FloatType fAxial, const Vector3D& relativeVelocityDirection, const Vector3D& normalCamber);

    FloatType computeBodyForcePerturbationScaling(size_t i, size_t j, size_t k, FloatType timePhysical) const;
    
};
