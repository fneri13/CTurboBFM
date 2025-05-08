#pragma once

#include "types.hpp"
#include "CMesh.hpp"
#include "Config.hpp"
#include "CFluid.hpp"
#include "commonFunctions.hpp"

// Base class handling BFM source terms computation
class CSourceBFMBase {
public:

    CSourceBFMBase(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : _config(config), _fluid(fluid), _mesh(mesh) {};

    virtual ~CSourceBFMBase() = default;  

    virtual StateVector computeSource(size_t i, size_t j, size_t k, const StateVector& primitive);

protected:
    const Config& _config;
    
    const CFluid& _fluid;
    
    const CMesh& _mesh;

    virtual StateVector computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive);

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive);

    Vector3D _point;
    FloatType _radius;
    FloatType _theta;
    FloatType _blockage;
    Vector3D _velocityCartesian;
    Vector3D _relativeVelocityCartesian;
    Vector3D _velocityCylindrical;
    FloatType _omega;
    Vector3D _dragVelocityCylindrical;
    Vector3D _relativeVelocityCylindric;
    FloatType _numberBlades;
    FloatType _pitch;
    FloatType _normalCamberAxial;
    FloatType _normalCamberRadial;
    FloatType _normalCamberTangential;
    Vector3D _normalCamberCylindric;
    FloatType _deviationAngle;
    Vector3D _inviscidForceDirectionCylindrical;
    Vector3D _viscousForceDirectionCylindrical;

    void computeFlowState(size_t i, size_t j, size_t k, const StateVector& primitive);

    FloatType computeDeviationAngle(Vector3D relativeVelocity, Vector3D normalCamber);

    Vector3D computeInviscidForceDirection(const Vector3D& relativeVelocity, const Vector3D& normalCamber);

    FloatType computeTangentialComponent(FloatType fAxial, const Vector3D& relativeVelocityDirection, const Vector3D& normalCamber);
    
};
