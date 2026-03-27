#pragma once

#include "types.hpp"
#include "mesh.hpp"
#include "config.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"
#include "math_utils.hpp"
#include "input_table.hpp"


class SourceBFMBase {
public:

    SourceBFMBase(const Config &config, const FluidBase &fluid, const Mesh &mesh);

    virtual ~SourceBFMBase() = default;  

    /** Compute the total source term: blockage + body force terms */
    virtual StateVector computeTotalSource(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce, 
        Matrix3D<Vector3D> &viscousForce, 
        Matrix3D<FloatType> &deviationAngle, 
        FloatType &timePhysical, 
        FlowSolution &conservativeVars,
        Matrix3D<FloatType> &timestep,
        FloatType &currentPhysicalTime);

protected:
    /** Thollet formulation */
    virtual StateVector computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive);

    virtual StateVector computeBodyForceSource(
        size_t i, 
        size_t j, 
        size_t k, 
        const StateVector& primitive, 
        Matrix3D<Vector3D> &inviscidForce, 
        Matrix3D<Vector3D> &viscousForce, 
        FlowSolution &conservativeVars,
        FloatType &dt,
        FloatType &timePhysical);
    
    void computeFlowState(
        size_t i, size_t j, size_t k, 
        const StateVector& primitive, 
        FlowSolution &conservativeVars, 
        FloatType &timePhysical);

    FloatType computeDeviationAngle(Vector3D relativeVelocity, Vector3D normalCamber);
    
    /** @brief The direction is obtained removing from the relative velocity the component 
     * aligned in the normal camber direction. 
     * Furthermore, any component out of the blade-to-blade plane is removed.
     */
    Vector3D computeInviscidForceDirection(
        const Vector3D& relativeVelocity, 
        const Vector3D& normalCamber, 
        size_t i, 
        size_t j, 
        size_t k);

    FloatType computeTangentialComponent(
        FloatType fAxial, 
        const Vector3D& relativeVelocityDirection, 
        const Vector3D& normalCamber);
    
    /** @brief Perturbation defined as forcePerturbed = forceBase * (1 + perturbationFactor)  */
    FloatType computePerturbationFactor(size_t i, size_t j, size_t k, FloatType timePhysical) const;


protected:
    const Config& _config;
    const FluidBase& _fluid;
    const Mesh& _mesh;
    
    bool _isBlockageActive {false};
    FloatType _radius;
    FloatType _theta;
    FloatType _blockage, _isBladePresent;
    FloatType _omega;
    FloatType _bladesNumber;
    FloatType _tangentialPitch, _staggeredPitch;
    FloatType _velMeridional;
    Vector3D _velCartesian, _velCylindrical;
    Vector3D _relVelCartesian, _dragVelCylindrical, _relVelCylindric;
    FloatType _normalCamberAxial, _normalCamberRadial, _normalCamberTangential;
    Vector3D _normalCamberCylindric;
    FloatType _deviationAngle;
    FloatType _gasPathAngle, _flowAngle, _leanAngle, _metalAngle;
    Vector3D _inviscidForceDirCylindrical, _inviscidForceDirCartesian;
    Vector3D _viscousForceDirCylindrical, _viscousForceDirCartesian;
    
    bool _isPerturbationActive {false};
    FloatType _perturbationRadius = 0.0;
    Vector3D _perturbationCenter = {0.0, 0.0, 0.0};
    FloatType _perturbationScaling = 1.0;
    FloatType _perturbationTimeStart = 0.0;
    FloatType _perturbationTimeDuration = 0.0;

    bool _isStalledBfmActive {false};
};
