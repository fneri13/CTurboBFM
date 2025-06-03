#include "CSourceBFMLiftDrag.hpp"

StateVector CSourceBFMLiftDrag::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    computeRelativeFlowAngle(i, j, k, primitive);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    // StateVector inviscidComponent = computeInviscidComponentGongStyle(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    return inviscidComponent + viscousComponent;
}

void CSourceBFMLiftDrag::computeRelativeFlowAngle(size_t i, size_t j, size_t k, const StateVector& primitive) {
    FloatType meridionalVelocityMag = std::sqrt(_relativeVelocityCylindric.x()*_relativeVelocityCylindric.x() + _relativeVelocityCylindric.y()*_relativeVelocityCylindric.y());

    // keep the convenetion used when computing beta_0. When the relative velocity is negative, the angle is negative
    _relativeFlowAngle = std::atan2(_relativeVelocityCylindric.z(), meridionalVelocityMag);
}


StateVector CSourceBFMLiftDrag::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType solidity = _mesh.getInputFields(FieldNames::SOLIDITY, i, j, k);
    FloatType hParameter = _mesh.getInputFields(FieldNames::LIFT_DRAG_H_PARAMETER, i, j, k);
    FloatType beta0 = _mesh.getInputFields(FieldNames::LIFT_DRAG_BETA_0, i, j, k);

    // compute the magnitude of the inviscid force. Positive when pushing, negative when pulling
    FloatType forceMag = 2.0 * M_PI * solidity * _relativeVelocityCylindric.dot(_relativeVelocityCylindric) / hParameter;    
    FloatType deltaBeta = std::abs(beta0 - _relativeFlowAngle);

    if (_deviationAngle < 0) {
        deltaBeta *= -1;
    }

    forceMag *= deltaBeta;
    
    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

StateVector CSourceBFMLiftDrag::computeInviscidComponentGongStyle(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType hParameter = _mesh.getInputFields(FieldNames::LIFT_DRAG_H_PARAMETER, i, j, k);
    FloatType knTurning = _mesh.getInputFields(FieldNames::LIFT_DRAG_KN_TURNING, i, j, k);

    // compute the magnitude of the inviscid force. Positive when pushing, negative when pulling
    FloatType forceMag = knTurning * _relativeVelocityCylindric.dot(_relativeVelocityCylindric) / hParameter * _deviationAngle;
    
    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

StateVector CSourceBFMLiftDrag::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    FloatType kp_etaMax = _mesh.getInputFields(FieldNames::LIFT_DRAG_KP_ETA_MAX, i, j, k);
    FloatType beta_etaMax = _mesh.getInputFields(FieldNames::LIFT_DRAG_BETA_ETA_MAX, i, j, k);
    FloatType solidity = _mesh.getInputFields(FieldNames::SOLIDITY, i, j, k);
    FloatType hParameter = _mesh.getInputFields(FieldNames::LIFT_DRAG_H_PARAMETER, i, j, k);
    
    FloatType kp = kp_etaMax + 2.0 * M_PI * solidity * std::pow(_relativeFlowAngle - beta_etaMax, 2);
        
    FloatType forceMag = _relativeVelocityCylindric.dot(_relativeVelocityCylindric) * kp / hParameter;

    Vector3D forceCylindrical = _viscousForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}